#!/usr/bin/env python

################################################################################
# SETUP
################################################################################
# Load required modules
import argparse
import sys
import os 
from collections import defaultdict
from itertools import product, combinations

import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score

# Load our modules
from munk.io import get_logger
from homolog_mapping_predictor import predict_sl_from_homolog_mapping as predict
from homolog_mapping_predictor import (
    load_homolog_mapping_file, 
    load_genetic_interactions_table, 
    compute_fpr)

################################################################################
# MAIN
################################################################################
# Parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-saf', '--species_A_files', type=str, required=True, nargs='*')
    parser.add_argument('-san', '--species_A_names', type=str, required=True, nargs='*') 
    parser.add_argument('-sap', '--species_A_ppis', type=str, required=True, nargs='*') 
    parser.add_argument('-sbf', '--species_B_files', type=str, required=True, nargs='*')
    parser.add_argument('-sbn', '--species_B_names', type=str, required=True, nargs='*') 
    parser.add_argument('-sbp', '--species_B_ppis', type=str, required=True, nargs='*') 
    parser.add_argument('-sn', '--species_names', type=str, required=True, nargs=2)
    parser.add_argument('-hmf', '--homolog_mapping_file', type=str, required=True)
    parser.add_argument('-of', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbose', type=int, required=False, default=1)
    parser.add_argument('-s', '--sample_negs', action='store_true')
    return parser.parse_args()

# Run
def run( args ):
    # Set up logging
    logger = get_logger(args.verbose)
    logger.info('* Loading input files...')
    
    # Load the homologs
    logger.info('\t- Loading homolog mapping...')
        
    geneAToHomologs, geneBToHomologs = \
        load_homolog_mapping_file(args.homolog_mapping_file)

    # Load the genetic interactions
    logger.info('\t- Loading genetic interaction files...')

    if args.sample_negs:
        logger.info('\t- Sampling Non-SLs...')
        
    As, Bs = [], []
    for species_A_file, species_A_ppi in \
        zip(args.species_A_files, args.species_A_ppis):
        As.append(load_genetic_interactions_table(
                    species_A_file, species_A_ppi, args.sample_negs))
    for species_B_file, species_B_ppi in \
        zip(args.species_B_files, args.species_B_ppis):
        Bs.append(load_genetic_interactions_table(
                    species_B_file, species_B_ppi, args.sample_negs))
    
    # Make predictions
    logger.info('* Making predictions for each pair of datasets...')
        
    items = []
    for (A_name, (A_GI, A_SL, A)), (B_name, (B_GI, B_SL, B)) in \
        product(zip(args.species_A_names, As), zip(args.species_B_names, Bs)):
        # Simple progress
        if args.verbose > 0:
            logger.info('\t- %s vs. %s' % (A_name, B_name))
            
        # Predict B from A
        B_edges, B_true, B_pred = predict((A_GI, A_SL, A),
                                          (B_GI, B_SL, B),
                                          geneBToHomologs,
                                          args.verbose)
        items.append({
            "Dataset A (Species)": '%s (%s)' % (A_name, args.species_names[0]),
            "Dataset B (Species)": '%s (%s)' % (B_name, args.species_names[1]),
            "Precision (True Positive Rate)": precision_score(B_true, B_pred),
            "Recall": recall_score(B_true, B_pred),
            "False Positive Rate": compute_fpr(B_true, B_pred),
            "F1 Score": f1_score(B_true, B_pred),
        })
        
        # Predict A from B
        A_edges, A_true, A_pred = predict((B_GI, B_SL, B),
                                          (A_GI, A_SL, A),
                                          geneAToHomologs,
                                          args.verbose)
        items.append({
            "Dataset A (Species)": '%s (%s)' % (B_name, args.species_names[1]),
            "Dataset B (Species)": '%s (%s)' % (A_name, args.species_names[0]),
            "Precision (True Positive Rate)": precision_score(A_true, A_pred),
            "Recall": recall_score(A_true, A_pred),
            "False Positive Rate": compute_fpr(A_true, A_pred),
            "F1 Score": f1_score(A_true, A_pred),
        })

    # Output to file
    logger.info('* Outputting to file...')
        
    df = pd.DataFrame(items)[ ['Dataset A (Species)', 'Dataset B (Species)', 
                               'Precision (True Positive Rate)', 'Recall', 
                               'False Positive Rate', 'F1 Score'] ]
    df.to_csv(args.output_file, index=0, sep='\t')

if __name__ == '__main__': run(parse_args())
