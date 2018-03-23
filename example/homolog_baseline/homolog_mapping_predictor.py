#!/usr/bin/env python

################################################################################
# SETUP
################################################################################
# Load required modules
import sys, os, argparse, pandas as pd
from collections import defaultdict
from itertools import product
from sklearn.metrics import precision_score, recall_score

# Load our modules
sys.path.append('../../src')
from i_o import get_logger
import random
import networkx as nx
from util import simple_two_core

################################################################################
# I/O AND PREDICTION
################################################################################
# Load the synthetic lethals
def load_genetic_interactions_table(gi_table_file, ppi_file,  sample_negs=False):
    logger = get_logger()
    GI = pd.read_csv(gi_table_file, sep='\t')

    # Explicitly filter out inconclusives...
    GI = GI[(GI['Category'] == 'SL') | (GI['Category'] == 'Non-SL')]
    SL = GI.loc[GI['Category'] == 'SL']
    pairs  = set( frozenset([u, v]) for u, v in zip(SL['Gene A'], SL['Gene B']))

    # If sample negatives...
    if sample_negs:
        assert(len(GI) == len(SL))

        nodes = set(simple_two_core(nx.read_edgelist(ppi_file)).nodes())
        non_sl_pairs = set()
        non_gi_data = []
        while len(non_sl_pairs) != len(pairs):
            p = frozenset(random.sample(nodes, 2))
            if p not in pairs and p not in non_sl_pairs:
                _p = tuple(p)
                non_sl_pairs.add( p )
                non_gi_data.append({'Category':'Non-SL', 'Score':1, 'Gene A': _p[0], 'Gene B': _p[1]})

        logger.info('\tAdding %s non-SL pairs for %s total interactions' % (len(non_sl_pairs), len(pairs)))
        pairs = pairs | set(non_sl_pairs)
        GI = pd.concat([GI, pd.DataFrame(non_gi_data)])
        SL = GI.loc[GI['Category'] == 'SL']

    return GI, SL, pairs

# Load a homolog mapping
def load_homolog_mapping_file(homolog_mapping_file):
    with open(homolog_mapping_file, 'r') as IN:
        homologs = set( frozenset(l.rstrip('\n').split()) for l in IN )
        geneAToHomologs = defaultdict(set)
        geneBToHomologs = defaultdict(set)
        for a_gene, b_gene in homologs:
            geneAToHomologs[a_gene].add( b_gene )

            geneBToHomologs[b_gene].add( a_gene )

    return geneAToHomologs, geneBToHomologs

def predict_sl_from_homolog_mapping(A_data, B_data, geneBToHomologs, verbose=1):
    A_GI, A_SL, A = A_data
    B_GI, B_SL, B = B_data

    # Make the predictions. We predict an SL between (u, v) in species B iff
    # (a) u and v have homologs (u', v') in species A; and, (b) (u', v') have an
    # SL.
    def hasSL(SLs, xs, ys):
        return any( frozenset([x, y]) in SLs for x, y in product(xs, ys) )

    B_true = [ int(c == 'SL') for c in B_GI['Category'] ]
    B_edges = [ (r['Gene A'], r['Gene B']) for i, r in B_GI.iterrows() ]
    B_pred = [ int(hasSL(A, geneBToHomologs[u], geneBToHomologs[v]))
               for u, v in B_edges ]

    return B_edges, B_true, B_pred

# Compute false positive rate for binary class predictions
def compute_fpr(true, pred):
    n_neg = len(true)-sum(true)
    return sum( 1. for x, y in zip(true, pred) if x == 0 and y == 1)/n_neg

################################################################################
# MAIN
################################################################################
# Parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sla', '--species_A_table_file', type=str, required=True)
    parser.add_argument('-slb', '--species_B_table_file', type=str, required=True)
    parser.add_argument('-hmf', '--homolog_mapping_file', type=str, required=True)
    parser.add_argument('-of', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbose', type=int, required=False, default=1)
    return parser.parse_args()

# Run
def run( args ):
    # Load the homologs
    geneAToHomologs, geneBToHomologs = load_homolog_mapping_file(args.homolog_mapping_file)
    
    # Load the genetic interactions
    A_GI, A_SL, A = load_genetic_interactions_table(args.species_A_table_file)
    B_GI, B_SL, B = load_genetic_interactions_table(args.species_B_table_file)

    logger = get_logger()
    if args.verbose > 0:
        logger.info('* Loaded genetic interactions...')
        logger.info('\t- %s interactions in species A (%s SLs)' % (len(A_GI), len(A)))
        logger.info('\t- %s interactions in species B (%s SLs)' % (len(B_GI), len(B)))

    # Make the predictions. We predict an SL between (u, v) in species B iff    # Load and predict
    B_edges, B_true, B_pred = predict_sl_from_homolog_mapping((A_GI, A_SL, A),
                                                              (B_GI, B_SL, B),
                                                              geneBToHomologs,
                                                              args.verbose)
    
    # Compute PR and AUC
    if args.verbose > 0:
        logger.info('* Evaluating results...') 
        logger.info('\tPrecision/TPR: %.3f' % precision_score(B_true, B_pred))
        logger.info('\tRecall: %.3f' % recall_score(B_true, B_pred))
        logger.info('\tFPR: %.3f' % compute_fpr(B_true, B_pred))

    ################################################################################
    # OUTPUT
    ################################################################################
    # Output predictions
    valToName = { 1: "SL", 0: "Non-SL" }
    items = [ {"Gene A": u, "Gene B": v, "Ground truth": valToName[gt], "Predicted": valToName[p]}
              for (u, v), gt, p in zip(B_edges, B_true, B_pred) ]
    df = pd.DataFrame(items)
    df.to_csv(args.output_file, sep='\t', index=False)

if __name__ == '__main__': run(parse_args())
