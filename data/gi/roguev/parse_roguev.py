#!/usr/bin/env python

# Load required modules
import sys, os, argparse
sys.path.append('..')
from gi_io import construct_gi_dataframe, output_gi_dataframe

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-o', '--output_prefix', type=str, required=True)
parser.add_argument('--sl_threshold', type=float, default=-2.5,
                    help='Default of -2.5 comes from the Data Quality section in the Supplement of Roguev, et al. (Science, 2008) supplementary information.')
parser.add_argument('--uncertainty_threshold', type=float, default=-1,
                    help='Default of -1 comes from the page 4 of the Collins, et al. (Nature, 2007) Supplement in the figure caption.')
parser.add_argument('-of', '--output_format', type=str, default='tsv', choices=['npy', 'tsv'])
args = parser.parse_args(sys.argv[1:])

# Load the E-MAP
pairToScore = dict()
n_discrepancies = 0
with open(args.input_file, 'r') as IN:
    arrs = [ l.rstrip('\n\r').split('\t') for l in IN ]
    gene_names = [ n.split('(')[0].replace('"', '') for n in arrs.pop(0)[1:] ]
    for arr in arrs:
        geneA = arr[0].split('(')[0].replace('"', '')
        for geneB, score in zip(gene_names, arr[1:]):
            if score == '': continue
            score = float(score)
            pair = frozenset([geneA, geneB])

            # If we measure the two genes twice and they're different,
            # we make the min by absolute value (and count how often this occurs)
            if pair in pairToScore:
                if pairToScore[pair] != score:
                    pairToScore[pair] = score if abs(score) < abs(pairToScore[pair]) else pairToScore[pair]
                    n_discrepancies += 1
            else:
                pairToScore[pair] = score

print '* Loaded %s interactions...' % len(pairToScore)
print '\t- Found %s discrepancies (pairs measured twice that were not the same)' % n_discrepancies
                
# Construct a dataframe and output to file
df = construct_gi_dataframe(pairToScore, args.sl_threshold, args.uncertainty_threshold)
output_gi_dataframe(df, args.output_prefix, args.output_format)
