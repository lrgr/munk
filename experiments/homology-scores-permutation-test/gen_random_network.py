import argparse
from collections import defaultdict
from time import time

import networkx as nx
import numpy as np
from sklearn.externals import joblib
from sklearn.externals.joblib import Parallel, delayed

from munk import util
from munk import regularized_laplacian, rkhs_factor
from munk.io import get_logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--edgelist', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-nt', '--n_tries', type=int, required=False, default=50)
    return parser.parse_args()

def shuffle(arr):
    ''' Copy and shuffle given list'''
    arr_copy = arr[:]
    np.random.shuffle(arr_copy)
    return arr_copy

def get_degree_to_names(tuples):
    '''
    Returns dictionary of degrees to node names from given list of
    (name, deg) tuples
    '''
    d = defaultdict(list)
    for  name, deg in tuples:
        d[deg].append(name)
    return d

def get_degree_preserving_relabelling(new_name_dict, unlabeled_G):
    '''
    Returns name mapping given between degree to node names dictionary
    and unlabeled graph G
    '''
    old_degree_dict = get_degree_to_names(unlabeled_G.degree)
    mapping_tuples = []
    for deg, new_names in new_name_dict.items():
        old_names = old_degree_dict[deg]
        mapping_tuples += (list(zip(old_names, shuffle(new_names))))
    return dict(mapping_tuples)

def perturbed_graph(seed_graph):
    '''
    Returns random graph with approximately the same degree
    distribution as given seed graph
    '''
    # Generate new graph with given degree sequence
    deg_sequence = shuffle(list(seed_graph.degree))
    new_G = nx.configuration_model(list(zip(*deg_sequence))[1])

    # Rename graph with nodes in seed_graph with corresponding degrees
    deg2names = get_degree_to_names(deg_sequence)
    mapping = get_degree_preserving_relabelling(deg2names, new_G)
    new_G = nx.relabel_nodes(new_G, mapping)

    # Remove self loops and parallel edges and return
    return util.simple_two_core(nx.Graph(new_G), verbose=False)

def safe_perturbed_graph(seed_graph, n_tries=10):
    if n_tries > 0:
        for i in range(n_tries):
            new_G = perturbed_graph(seed_graph)
            if len(not_shared_nodes(new_G, seed_graph)) == 0:
                return new_G, i + 1
    raise Exception()

def not_shared_nodes(G1, G2):
    return set(G1.nodes) ^ set(G2.nodes)

def gen_random_network_data(seed_graph, n_tries=10):
    log = get_logger()
    random_graph, tries = safe_perturbed_graph(seed_graph, n_tries)
    nodes = sorted(random_graph.nodes())
    D = regularized_laplacian(random_graph, nodes, lam=0.05)
    C = rkhs_factor(D)
    return dict(G=random_graph, D=D, C=C, nodes=nodes), tries

def main(args):
    log = get_logger()
    log.info('Loading edgelist from: %s', args.edgelist)
    seed_network = \
        util.simple_two_core(nx.read_edgelist(args.edgelist, encoding='ascii'),
                             verbose=False)

    log.info('Generating network...')
    t_start = time()
    data, tries = gen_random_network_data(seed_network, args.n_tries)
    t_end = time()
    elapsed = t_end - t_start
    log.info('Random network generated in %.2f, in %d tries', elapsed, tries)

    log.info('Saving random network to: %s', args.output)
    joblib.dump(data, args.output)


if __name__ == '__main__':
    main(parse_args())
