#!/usr/bin/env

# Load required modules
import sys, os, argparse
sys.path.append('..')
from gi_io import construct_gi_dataframe, output_gi_dataframe

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-o', '--output_prefix', type=str, required=True)
parser.add_argument('--accepted_perturbations', type=str, required=False, nargs='*',
                    default=['deletion'], choices=['deletion', 'damp'])
parser.add_argument('--sl_threshold', type=float, default=-3,
                    help='Default of -3 comes from figure caption on page 4 of Collins, et al. (Nature, 2007) supplementary information.')
parser.add_argument('--uncertainty_threshold', type=float, default=-1,
                    help='Default of -1 comes from the page 4 of the Collins, et al. (Nature, 2007) Supplement in the figure caption.')
parser.add_argument('-of', '--output_format', type=str, default='tsv', choices=['npy', 'tsv'])
args = parser.parse_args(sys.argv[1:])

# Load the raw EMAP file
accepted_perturbations = set([ p.upper() for p in args.accepted_perturbations ])
pairToScore = dict()
with open(args.input_file, 'r') as IN:
    arrs = [ l.rstrip('\n').split('\t') for l in IN ]
    gene_names = arrs.pop(0)[3:]
    perturbations = arrs.pop(0)[3:]
    arrs.pop(0) # skip the third line
    for arr in arrs:
        geneA, perturbation = arr[0], arr[1]
        if perturbation in accepted_perturbations:
            for i, (geneB, score) in enumerate(zip(gene_names, arr[3:])):
                if perturbations[i] in accepted_perturbations and score != 'NaN':
                    # if there is a duplicate pair, it takes the most recent (but
                    # they should be identical)
                    score = float(score)
                    if frozenset([geneA, geneB]) in pairToScore:
                        assert(pairToScore[frozenset([geneA, geneB])] == score)
                    pairToScore[frozenset([geneA, geneB])] = score
                    
# Construct a dataframe and output to file
df = construct_gi_dataframe(pairToScore, args.sl_threshold, args.uncertainty_threshold)
output_gi_dataframe(df, args.output_prefix, args.output_format)
