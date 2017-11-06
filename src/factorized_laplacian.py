
# Load required modules

import argparse
import os
import sys
import time

import numpy as np
import networkx as nx
import scipy
from sklearn.externals import joblib

import logging_utils
import util
from network import regularized_laplacian, rkhs_factor

##########################################################################
# MAIN
##########################################################################
# Command-line argument parser


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--edge_file', type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-p', '--dimension_percentage', type=float, default=0.1)
    parser.add_argument('-l', '--lam', type=float, default=0.05)
    parser.add_argument('-df', '--diffusion_file', type=str)
    return parser


# for use in development 
def setup_args(ef='data/collins/collins-confident-ppi-edges.txt',
               of='output/findr_features/collins-ppi-findr.pkl'):
    args = argparse.Namespace()
    args.edge_file = ef
    args.output_file = of 
    args.dimension_percentage = 0.1
    args.lam = 0.05
    return args

def run(args):
    logger = logging_utils.getLogger()

    # Sanity Checks
    dimension_percentage = args.dimension_percentage
    assert 0.0 <= dimension_percentage and dimension_percentage <= 1

    G = nx.read_edgelist(args.edge_file, encoding='ascii')
    G = util.simple_two_core(G)

    nodes = sorted(G.nodes())

    L = np.array(nx.laplacian_matrix(G, nodelist=nodes).todense())

    # compute regularized laplacian
    logger.info('Computing regularized Laplacian')
    D = regularized_laplacian(L, args.lam)

    if args.diffusion_file:
        output = dict(D = D, nodes = nodes)
        joblib.dump(output, args.diffusion_file)

    # extract RKHS vectors
    logger.info('Extracting RKHS features')
    findr = rkhs_factor(D)

    reduced_dims = int(len(nodes) * dimension_percentage)

    # Output to file
    output = dict(params=vars(args), G=G, D=findr,
                  W=findr[:,-reduced_dims:], X=findr[:,-reduced_dims:], nodes=nodes)
    joblib.dump(output, args.output_file)


if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
