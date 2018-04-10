#!/usr/bin/env python

################################################################################
# SETUP
################################################################################
# Load required modules
import sys, os, argparse, logging, pandas as pd, numpy as np
from sklearn.externals import joblib
from collections import Counter

# Load our modules
from constants import *
from handl.io import get_logger

################################################################################
# MAIN
################################################################################
# Parse command-line arguments
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ef', '--embedding_file', type=str, required=True)
    parser.add_argument('-gif', '--gi_file', type=str, required=True)
    parser.add_argument('-of', '--output_file', type=str, required=True)
    parser.add_argument('-ff', '--feature_function', type=str, required=False,
                        default=ADD_FEATURES, choices=FEATURE_FUNCTIONS)
    parser.add_argument('-v', '--verbosity', type=int, required=False,
                        default=logging.INFO)
    parser.add_argument('--sinatra-featurize', action='store_true', default=False)
    parser.add_argument('--remove-landmarks', action='store_true', default=False)
    return parser

# Main
def run( args ):
    # Create logger
    logger = get_logger(args.verbosity)

    # Load the SLs
    logger.info('[Loading genetic interactions file]')
    df = pd.read_csv(args.gi_file, sep='\t', dtype={'Gene A': str, 'Gene B': str, 'Score': float})
    pairToOutcome = dict( (frozenset([u, v]), c) for u, v, c in zip(df['Gene A'], df['Gene B'], df['Category']) if c != INCONCLUSIVE )

    if args.sinatra_featurize:
        # Load only SLs for SINaTRA features
        pairs = sorted( p for p, c in pairToOutcome.items() if c == SL )
        logger.info('- SINaTRA features: loaded %s SLs' % len(pairs))
    else:
        # We include SLs and non-SLs
        pairs = sorted( p for p, c in pairToOutcome.items() )
        
        # Do some light reporting
        class_counter = Counter(pairToOutcome.values())
        n_sl = class_counter[SL]
        n_non_sl = class_counter[NON_SL]
        logger.info('- Loaded %s interactions (%s SL and %s non-SL)' % (len(pairs), n_sl, n_non_sl))

    # Load the embedding
    logger.info('[Loading HANDL embedding]')
    obj = joblib.load(args.embedding_file)
    embedding = obj.get('X')
    nodes = obj.get('nodes')
    landmark_nodes = set(obj.get('landmarks'))
    node_set = set(nodes)
    n_nodes = len(nodes)
    nodeToIndex = dict(zip(nodes, range(n_nodes)))
    logger.info('- %s nodes' % n_nodes)
    logger.info('- %s landmarks' % len(landmark_nodes))
        
    # Restrict to pairs also in the network
    pairs = [ p for p in pairs if all([ u in node_set for u in p ]) and len(p) == 2 ]
    logger.info('- Restricting to %s pairs in the network' % len(pairs))
    
    if args.remove_landmarks:
        pairs = [ p for p in pairs if all([ u not in landmark_nodes for u in p ]) ]
        logger.info('- Restricting to %s pairs in the network without landmarks' % len(pairs))
    
    if args.sinatra_featurize:
        # For SINaTRA features, we randomly sample a set of
        # non-SLs from nodes in the network
        import random
        sl_pairs = list(pairs)
        non_sl_pairs = set()
        while len(non_sl_pairs) != len(sl_pairs):
            if args.remove_landmarks:
                p = frozenset(random.sample(node_set-landmark_nodes, 2))
            else:
                p = frozenset(random.sample(nodes, 2))
            if p not in sl_pairs and p not in non_sl_pairs:
                non_sl_pairs.add( p )
                
        pairToOutcome.update( (p, NON_SL) for p in non_sl_pairs )
        pairs = sl_pairs + list(non_sl_pairs)

        logger.info('\tAdding %s non-SL pairs for %s total interactions' % (len(non_sl_pairs), len(pairs)))

    # Determine the feature function
    if args.feature_function == ADD_FEATURES:
        def feature_function(p):
            u, v = sorted(p)
            return embedding[nodeToIndex[u]] + embedding[nodeToIndex[v]]
        
    elif args.feature_function == MEAN_FEATURES:
        def feature_function(p):
            u, v = sorted(p)
            return (embedding[nodeToIndex[u]] + embedding[nodeToIndex[v]])/2.
    else:
        raise NotImplementedError('Feature function "%s" not implemented' % args.feature_function)
        
    # Construct the features and outcomes and output to file
    n_pairs = len(pairs)
    X = np.array([ feature_function(p) for p in pairs ])
    y = np.array([ pairToOutcome[p] == SL for p in pairs ])

    # Output to file
    output = dict(X=X, y=y, pairs=pairs, params=vars(args))
    joblib.dump(output, args.output_file)
    
if __name__ == '__main__':
    run( get_parser().parse_args(sys.argv[1:]) )
