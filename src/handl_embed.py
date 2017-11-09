#!/usr/bin/env Python 

import argparse
import os
import sys
import time

import networkx as nx
import numpy as np
import scipy
from sklearn.externals import joblib

import logging_utils

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--source_rkhs_file', type=str, required=True)
    parser.add_argument('-t', '--target_laplacian_file', type=str, required=True)
    parser.add_argument('-hf', '--homologs_file', type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-n', '--n_landmarks', type=int, required=False, default=400)
    return parser

def run(args):
    logger = logging_utils.getLogger()

    # read in source and target matrices
    data = joblib.load(args.source_rkhs_file)
    source_C = data['X']
    source_nodes = data['nodes']
    
    data = joblib.load(args.target_laplacian_file)
    target_nodes = data['nodes']
    target_D = data['D']

    # read in homologs
    # assume the order is source node - target node
    source_homs = []
    target_homs = []
    with open(args.homologs_file, 'r') as hom_f:
        for line in hom_f:
            source_hom, target_hom = line.split()
            if (source_hom in source_nodes) and (target_hom in target_nodes):
                source_homs.append(source_hom)
                target_homs.append(target_hom)

    n_landmarks = args.n_landmarks

    # Sanity check for number of landmarks
    assert(n_landmarks <= len(source_homs))
   
    source_homolog_indices = [source_nodes.index(node) for node in source_homs]
    target_homolog_indices = [target_nodes.index(node) for node in target_homs]
    source_landmarks = source_homs[:n_landmarks]
    target_landmarks = target_homs[:n_landmarks]

    source_landmark_indices = source_homolog_indices[:n_landmarks]
    source_non_landmark_indices = source_homolog_indices[n_landmarks:]
    target_landmark_indices = target_homolog_indices[:n_landmarks]
    target_non_landmark_indices = target_homolog_indices[n_landmarks:]

    # select corresponding rows in each vector set and regress
    target_handl_C = np.linalg.pinv(source_C[source_landmark_indices,:]).dot(
        target_D[target_landmark_indices,:])

    data = dict(target_handl_C = target_handl_C.T, 
                source_C = source_C,
                D = source_C.dot(target_handl_C), 
                target_nodes = target_nodes, 
                source_nodes = source_nodes,
                target_landmarks = target_landmarks,
                target_landmark_indices = target_landmark_indices,
                target_non_landmark_indices = target_non_landmark_indices,
                source_landmarks = source_landmarks,
                source_landmark_indices=source_landmark_indices,
                source_non_landmark_indices = source_non_landmark_indices)
    joblib.dump(data, args.output_file)


    logger.info('Num landmarks used for embedding: %d', n_landmarks)
    logger.info('Source embedding  shape: %s', str(source_C.shape))
    logger.info('Pre-embedding target shape: %s', str(target_D.shape))
    logger.info('Post-embedding targets hape: %s', str(target_handl_C.T.shape))


if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
