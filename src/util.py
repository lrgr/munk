#!/usr/bin/env python

# Load required modules
import numpy as np
import networkx as nx 
import pandas as pd

from logging_utils import getLogger

# Vector operations

def simplify_graph(G):
    ''' 
    Returns the simple/strict graph corresponding to given graph
    (Removes self loops from G and returns largest connected component)
    '''
    logger = getLogger()
    if (not nx.is_connected(G)):
        cc_list = list(nx.connected_component_subgraphs(G))
        cc_sizes = [len(x) for x in cc_list]
        largest_cc = max(cc_sizes)
        cc_sizes.remove(largest_cc)

        logger.warning('Network has %d connected components', len(cc_list))
        logger.warning('\tLargest is size %d and all the rest are %d or smaller', 
            largest_cc, max(cc_sizes))
        logger.warning('\tUsing largest connected component')

        G = max(cc_list, key=len)
    G.remove_edges_from(list(G.selfloop_edges()))
    return G


def simple_two_core(G):
    ''' Returns simplified, 2 cored of given graph '''

    logger = getLogger()

    # Get simple graph
    G = simplify_graph(G)

    # Compute 2 core
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    logger.info('PPI info - # Nodes: %d, # Edges: %d', num_nodes, num_edges)

    logger.info('Computing 2 core')
    G = nx.k_core(G, 2)
    num_2core_nodes = G.number_of_nodes()
    num_2core_edges = G.number_of_edges()
    logger.info('2 core info - # Nodes: %d, # Edges: %d', 
        num_2core_nodes, num_2core_edges)
    logger.info('2 core removed %d nodes and %d edges', 
        num_nodes - num_2core_nodes, num_edges - num_2core_edges)
    return G

def binarized_graphs_from_tsv(fp):
    # load tsv
    df = pd.DataFrame.from_csv(fp, sep='\t', header=0, index_col=None)
    
    sls_df = df[df['Category'] == 'SL'].values
    nonsls_df = df[df['Category'] == 'Non-SL'].values

    SLs = nx.Graph() 
    non_SLs = nx.Graph()


    SLs.add_weighted_edges_from(sls_df[:, :3])
    non_SLs.add_weighted_edges_from(nonsls_df[:, :3])

    return SLs, non_SLs

def graph_from_tsv(fp):
    # load tsv
    df = pd.DataFrame.from_csv(fp, sep='\t', header=0, index_col=None)
    
    G = nx.Graph()
    G.add_weighted_edges_from(df.values[:, :3])
    return G
