# handl.py Implementation of HANDL
# Load required modules
import argparse
import json
import random
import time

import networkx as nx
import numpy as np
import scipy as sp
from sklearn.externals import joblib

import util
from i_o import get_logger

###############################################################################
# utility to seperate scores
###############################################################################

###############################################################################
# DIFFUSION MEASURES
###############################################################################

def regularized_laplacian(G, nodes, lam):
    ''' 
    Computes regularized laplacian from networkx graph corresponding to
    given node ordering

    Note: index of regularized laplacian rows/columns correspond to indices of
    the sorted list of nodes in the given graph.
    '''
    L = np.array(nx.laplacian_matrix(G, nodelist=nodes).todense())
    return np.linalg.inv(np.eye( *np.shape(L) ) + (lam * L))

def rkhs_factor(D):
    ''' Computes RKHS embeddings for given kernel matrix '''
    e, v = sp.linalg.eigh(D)
    return v.dot(np.diag(np.sqrt(e)))

###############################################################################
# HANDL embedding
###############################################################################

def handl_embed_matrices(source_C, target_D, landmark_idxs):
    ''' 
    Computes HANDL embeddings of source and target matrices given corresponding
    indices of landmarks

    :param source_C: 2D array of rkhs vectors from source species
    :param source_D: 2D array of diffusion values from target species
    :param landmark_idxs: list of landmark tuples like (source_idx, target_idx)
j   
    :return: tuple of HANDL embeddings for source and target species
    '''
    source_idxs, target_idxs = zip(*landmark_idxs)
    target_C_hat = np.linalg.pinv(source_C[source_idxs,:]) \
                    .dot(target_D[target_idxs,:]).T

    return source_C, target_C_hat

def handl_embed_graphs(source_G, target_G, homologs, n_landmarks, lam=0.05, return_idxs=False):
    '''
    Computes HANDL embeddings of given source and target graphs with given
    list homologs and number of landmarks

    :param source_G: networkx graph corresponding to source species
    :param target_G: networkx graph corresponding to target species
    :param homologs: list of tuples of homologs like (source_node, target_node)
    :param n_landmarks: The number of landmarks to use for HANDL embedding
    :param lam: lambda value for regularized laplacian
    
    :return: HANDL embeddings for source and target species
    '''
    assert(n_landmarks <= len(homologs))

    source_nodes = util.sorted_nodes(source_G)
    target_nodes = util.sorted_nodes(target_G)

    t_start = time.time()
    source_D = regularized_laplacian(source_G, source_nodes, lam)
    t_end = time.time()
    source_laplacian_time = t_end - t_start

    t_start = time.time()
    source_C = rkhs_factor(source_D)
    t_end = time.time()
    source_rkhs_time = t_end - t_start 

    t_start = time.time()
    target_D = regularized_laplacian(target_G, target_nodes, lam)
    t_end = time.time()
    target_laplacian_time = t_end - t_start

    source_node2index = util.node_to_index(source_G)
    target_node2index = util.node_to_index(target_G)

    landmark_homs = homologs[:n_landmarks]
    landmark_idxs = [(source_node2index[s_n], target_node2index[t_n]) for 
                     s_n, t_n in landmark_homs]

    t_start = time.time()
    source_C, target_C_hat = handl_embed_matrices(source_C, target_D, landmark_idxs)
    t_end = time.time()
    handl_time = t_end - t_start

    runtimes = dict(source_regularized_laplacian_runtime=source_laplacian_time,
                    source_rkhs_factorization_runtime=source_rkhs_time,
                    target_regularized_laplacian_runtime=target_laplacian_time,
                    handl_embed_runtime=handl_time)

    if return_idxs:
        homolog_idxs = [(source_node2index[s_n], target_node2index[t_n]) for 
                        s_n, t_n in homologs] 
        return ((source_C, source_nodes), \
                (target_C_hat, target_nodes), \
                landmark_idxs, \
                homolog_idxs)
    else:
        return ((source_C, source_nodes), \
                (target_C_hat, target_nodes), \
                landmark_homs, runtimes)


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

###############################################################################
# MAIN
###############################################################################

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-se', '--source_edgelist', type=str, required=True)
    parser.add_argument('-te', '--target_edgelist', type=str, required=True)
    parser.add_argument('-hf', '--homolog_list', type=str, required=True)
    parser.add_argument('-so', '--source_output_file', type=str, required=True)
    parser.add_argument('-to', '--target_output_file', type=str, required=True)
    parser.add_argument('-sim', '--sim_scores_output_file', type=str, required=True)
    parser.add_argument('-lo', '--landmarks_output_file', type=str, required=True)
    parser.add_argument('-r', '--runtimes_file', type=str, required=True)
    parser.add_argument('-n', '--n_landmarks', type=int, required=False, default=400)
    return parser.parse_args()

def main(args):
    '''
    Compute HANDL embeddings given target and species edgelists and list of
    homologs.
    '''
    log = get_logger()

    log.info('Loading homologs list from %s', args.homolog_list)
    raw_homologs = util.read_homolog_list(args.homolog_list)

    log.info('Loading source edgelist from %s', args.source_edgelist)
    source_G = nx.read_edgelist(args.source_edgelist, encoding='ascii')
    log.info('Loading target edgelist from %s', args.target_edgelist)
    target_G = nx.read_edgelist(args.target_edgelist, encoding='ascii')

    log.info('Computing HANDL embeddings with %d landmarks', args.n_landmarks)
    n_landmarks = args.n_landmarks

    source_G = util.simple_two_core(source_G)
    target_G = util.simple_two_core(target_G)
    homologs = homologs_in_graphs(source_G, target_G, raw_homologs)

    t_start = time.time()
    source_data, target_data, landmarks, runtimes = \
        handl_embed_graphs(source_G, target_G, homologs, n_landmarks)
    t_end = time.time()
    total_time = t_end - t_start
    source_X, source_nodes = source_data
    target_X, target_nodes = target_data

    log.info('Saving source species embeddings to %s', args.source_output_file)
    util.save_embeddings(source_X, source_nodes, args.source_output_file)
    log.info('Saving target species embeddings to %s', args.target_output_file)
    util.save_embeddings(target_X, target_nodes, args.target_output_file)

    log.info('Source data shape {}'.format(source_X.shape))
    log.info('Target data shape {}'.format(target_X.shape))
    
    log.info('Saving HANDL similarity scores to %s', args.sim_scores_output_file)
    
    sim_scores_matrix = np.dot(source_X, target_X.T)
    joblib.dump(dict(X=sim_scores_matrix, A_nodes=source_nodes, 
                     B_nodes=target_nodes),
                args.sim_scores_output_file)

    log.info('Saving landmark list to %s', args.landmarks_output_file)

    with open(args.landmarks_output_file, 'w') as OUT:
        for a, b in landmarks:
            OUT.write('%s\t%s\n' % (a, b)) 

    log.info('Saving runtimes to %s', args.runtimes_file)
    with open(args.runtimes_file, 'w') as OUT:
        json.dump(dict(n_source_nodes = source_G.number_of_nodes(),
                       n_target_nodes = target_G.number_of_nodes(),
                       runtimes=runtimes,
                       total_time=total_time), OUT, indent=2)

    log.info('HANDL embedding complete!')

if __name__ == '__main__':
    main(parse_args())
