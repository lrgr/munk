
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
    parser.add_argument('-l', '--lam', type=float, default=0.05)
    parser.add_argument('-df', '--diffusion_file', type=str, required=True)
    return parser

def run(args):
    logger = logging_utils.getLogger()

    G = nx.read_edgelist(args.edge_file, encoding='ascii')
    G = util.simple_two_core(G)

    nodes = sorted(G.nodes())

    L = np.array(nx.laplacian_matrix(G, nodelist=nodes).todense())

    # compute regularized laplacian
    logger.info('Computing regularized Laplacian')
    D = regularized_laplacian(L, args.lam)

    output = dict(params=vars(args), D = D, nodes = nodes)
    joblib.dump(output, args.diffusion_file)

    # extract RKHS vectors
    logger.info('Extracting RKHS features')
    rkhs = rkhs_factor(D)

    # Output to file
    output = dict(params=vars(args), X = rkhs, nodes = nodes)
    joblib.dump(output, args.output_file)

if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
