import argparse

import networkx as nx
import numpy as np

import util

from collections import defaultdict
from sklearn.externals import joblib
from sklearn.externals.joblib import Parallel, delayed
from i_o import get_logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species_name', type=str, required=True)
    parser.add_argument('-n', '--n_graphs', type=int, required=True)
    parser.add_argument('-e', '--edgelist', type=str, required=True)
    parser.add_argument('-d', '--output_dir', type=str, required=True)
    parser.add_argument('-j', '--n_jobs', type=int, required=False, default=1)
    return parser.parse_args()

def shuffle(arr):
    ''' Copy and shuffle given list'''
    arr_copy = arr[:]
    np.random.shuffle(arr_copy)
    return arr_copy

def get_degree_to_names(tuples):
    ''' Returns dictionary of degrees to node names from given list of (name, deg) tuples '''
    d = defaultdict(list)
    for  name, deg in tuples:
        d[deg].append(name)
    return d

def get_degree_preserving_relabelling(new_name_dict, unlabeled_G):
    ''' Returns name mapping given between degree to node names dictionary and unlabeled graph G'''
    old_degree_dict = get_degree_to_names(unlabeled_G.degree)
    mapping_tuples = []
    for deg, new_names in new_name_dict.items():
        old_names = old_degree_dict[deg]
        mapping_tuples += (list(zip(old_names, shuffle(new_names))))
    return dict(mapping_tuples)

def perturbed_graph(seed_graph):
    ''' Returns random graph with approximately the same degree distribution as given seed graph'''
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
    for i in range(n_tries):
        new_G = perturbed_graph(seed_graph)
        if len(not_shared_nodes(new_G, seed_graph)) == 0:
            return new_G
    raise Exception()

def not_shared_nodes(G1, G2):
    return set(G1.nodes) ^ set(G2.nodes)

def gen_and_save_random_network(fp, seed_graph, n_tries=10):
	log = get_logger()
	random_graph = safe_perturbed_graph(seed_graph, n_tries)
	log.info('Saving random network to: %s', fp)
	joblib.dump(random_graph, fp)


def main(args):
	log = get_logger()
	log.info('Loading edgelist from: %s', args.edgelist)
	log.info('Saving random networks to directory: %s', args.output_dir)
	log.info('Generating %d networks', args.n_graphs)

	filenames = ['{}{}.pkl'.format(args.species_name, i) for i in range(args.n_graphs)]
	filepaths = ['{}/{}'.format(args.output_dir, filename) for filename in filenames]

	seed_network = util.simple_two_core(nx.read_edgelist(args.edgelist, encoding='ascii'), verbose=False)
	Parallel(n_jobs=args.n_jobs)(
			delayed(gen_and_save_random_network)(fp, seed_network) for fp in filepaths)

	with open('{}/{}-networks.tsv'.format(args.output_dir, args.species_name),'w') as OUT:
		for fname in filenames:
  			OUT.write("%s\n" % fname)

if __name__ == '__main__':
	main(parse_args())
