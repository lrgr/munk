#!/usr/bin/env Python

# Load required modules

import argparse
import os
import sys
import time

import networkx as nx
import numpy as np
import scipy
from sklearn.externals import joblib

import logging_utils

def announce(message):
    print time.strftime('%H:%M:%S'),message
    sys.stdout.flush()

##########################################################################
# MAIN
##########################################################################
# Command-line argument parser


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sv', '--static_vector_file', type=str, required=True)
    parser.add_argument('-df', '--diffusion_file', type=str, required=True)
    parser.add_argument('-hf', '--homologs_file', type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-n', '--n_landmarks', type=int, required=False, default=400)
    return parser

# for use in development 
def setup_args(sv='output/features/biogrid-sc-ppi-findr.pkl',
               df='output/features/biogrid-sp-diffusion.pkl',
               h = '../../data/homologs/sc-sp/sc-sp-homologs.txt',
               o = 'output/features/biogrid-sp-ppi-aligned.pkl'
               ):
    args = argparse.Namespace()
    args.static_vector_file = sv
    args.diffusion_file = df
    args.homologs_file = h
    args.output_file = o
    return args

def run(args):
    logger = logging_utils.getLogger()

    # read in static vector file and its nodelist
    data = joblib.load(args.static_vector_file)
    static_X = data['X']
    static_nodes = data['nodes']
    
    data = joblib.load(args.diffusion_file)
    aligned_nodes = data['nodes']
    alignment_metric = data['D']

    # read in homologs
    # assume the order is static node - aligned node
    # this means sc - sp
    static_homs = []
    aligned_homs = []
    with open(args.homologs_file, 'r') as hom_f:
        line = hom_f.readline()
        while line:
            static_hom, aligned_hom = line.split()
            if (static_hom in static_nodes) and (aligned_hom in aligned_nodes):
                static_homs.append(static_hom)
                aligned_homs.append(aligned_hom)
            line = hom_f.readline()

    n_homologs = args.n_landmarks

    # because we are using only a fraction of the findr feature vectors
    # we want to make sure we havent used more landmarks than feature dimension
    assert(n_homologs < static_X.shape[1])
   
    static_indices = [static_nodes.index(node) for node in static_homs]
    aligned_indices = [aligned_nodes.index(node) for node in aligned_homs]
    static_landmark_homs = static_homs[:n_homologs]
    aligned_landmark_homs = aligned_homs[:n_homologs]

    static_landmark_indices = static_indices[:n_homologs]
    static_non_landmark_hom_indices = static_indices[n_homologs:]
    aligned_landmark_indices = aligned_indices[:n_homologs]
    aligned_non_landmark_hom_indices = aligned_indices[n_homologs:]

    # select corresponding rows in each vector set and regress
    aligned_X = np.linalg.pinv(static_X[static_landmark_indices,:]).dot(
        alignment_metric[aligned_landmark_indices,:])

    # use aligned_X.T
    data = dict(aligned_X = aligned_X.T, 
                static_X = static_X,
                D = static_X.dot(aligned_X), 
                nodes = aligned_nodes, 
                static_nodes = static_nodes,
                aligned_landmark_homs = aligned_landmark_homs,
                aligned_landmark_indices=aligned_landmark_indices,
                aligned_non_landmark_homs_indices = aligned_non_landmark_hom_indices,
                static_landmark_homs = static_landmark_homs,
                static_landmark_indices=static_landmark_indices,
                static_non_landmark_homs_indices = static_non_landmark_hom_indices)
    joblib.dump(data, args.output_file)


    logger.info('Num homologs used for alignment: %d', n_homologs)
    logger.info('Pre-alignment shape: %s', str(static_X.shape))
    logger.info('Post-alignment shape: %s', str(aligned_X.T.shape))
    logger.info('Num aligned nodes: %d', len(aligned_nodes))
    logger.info('alignment metric shape: %s', str(alignment_metric.shape))



if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
