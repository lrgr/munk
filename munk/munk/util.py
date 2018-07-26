#!/usr/bin/env python

# Load required modules
import sys, os, numpy as np

import networkx as nx
from sklearn import metrics
from sklearn.externals import joblib

from munk.io import get_logger

def simplify_graph(G, verbose=True):
    '''
    Returns the simple/strict graph corresponding to given graph
    (Removes self loops from G and returns largest connected component)
    '''
    logger = get_logger()
    if (not nx.is_connected(G)):
        cc_list = list(nx.connected_component_subgraphs(G))
        cc_sizes = [len(x) for x in cc_list]
        largest_cc = max(cc_sizes)
        cc_sizes.remove(largest_cc)
        if verbose:
            logger.warning('Network has %d connected components', len(cc_list))
            logger.warning('\tLargest is size %d and all the rest are %d or smaller',
                largest_cc, max(cc_sizes))
            logger.warning('\tUsing largest connected component')

        G = max(cc_list, key=len)
    G.remove_edges_from(G.selfloop_edges())
    return G

def simple_two_core(G, verbose=True):
    ''' Returns simple, 2 core of given graph '''

    logger = get_logger()
    # Get simple graph
    G = simplify_graph(G, verbose)
    # Compute 2 core
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    if verbose:
        logger.info('PPI info - # Nodes: %d, # Edges: %d', num_nodes, num_edges)
        logger.info('Computing 2 core')
    G = nx.k_core(G, 2)
    num_2core_nodes = G.number_of_nodes()
    num_2core_edges = G.number_of_edges()
    if verbose:
        logger.info('2 core info - # Nodes: %d, # Edges: %d',
            num_2core_nodes, num_2core_edges)
        logger.info('2 core removed %d nodes and %d edges',
            num_nodes - num_2core_nodes, num_edges - num_2core_edges)
    return G

def read_homolog_list(fp):
    ''' Returns tuple of homologs from two column homolgs tsv file '''
    with open(fp, 'r') as f:
        return [tuple(l.split()) for l in f]

def homologs_in_graphs(G1, G2, homologs):
    '''
    Computes list of homologs from a pair of graphs and a given list of
    homolog candidates
    '''
    G1_nodes = set(G1.nodes())
    G2_nodes = set(G2.nodes())
    valid_homologs = [(h1, h2) for h1, h2 in homologs \
                        if h1 in G1_nodes and h2 in G2_nodes]
    return valid_homologs

# Functions to easily convert between nodes and indices
def node_to_index(obj, dtype='graph'):
    if dtype == 'graph':
        return dict((n,i) for i, n in enumerate(sorted_nodes(obj)))
    elif dtype == 'list':
        return dict((n,i) for i, n in enumerate(obj))
    else:
        raise NotImplementedError()

def index_to_node(obj, dtype='graph'):
    if dtype == 'graph':
        return dict(enumerate(sorted_nodes(obj)))
    elif dtype == 'list':
        return dict(enumerate(obj))

def sorted_nodes(G):
    return sorted(G.nodes())
