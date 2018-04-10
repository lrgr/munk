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

import handl.util as util
from handl.io import get_logger

###############################################################################
# utility to seperate scores
###############################################################################

def separate_scores(scores, landmark_pair_idxs, homolog_pair_idxs):
    ''' 
    Separate scores into the following 5 categories:
        1. Landmark - Landmark pair scores
        2. Off diagonal entries in the landmark-landmark submatrix
        3. Landmark- (non-landmark) entries in rows and columns that correspond
           to one landmark
        4. Homolog-homolog pairs
        5. Other pairs
    '''
    source_landmark_idxs, target_landmark_idxs = zip(*landmark_pair_idxs)
    source_homolog_idxs, target_homolog_idxs = zip(*homolog_pair_idxs)
    landmark_mask = np.zeros_like(scores, dtype=bool)
    
    source_landmark_target_all_mask = np.zeros_like(scores, dtype=bool)
    source_landmark_target_all_mask[source_landmark_idxs, :]  = True
    source_all_target_landmark_mask = np.zeros_like(scores, dtype=bool)
    source_all_target_landmark_mask[:, target_landmark_idxs]  = True
    landmark_landmark_mask = np.zeros_like(scores, dtype=bool) 
    landmark_landmark_mask[source_landmark_idxs, target_landmark_idxs] = True
    
    # Obtain landmark-landmark pairs
    L_L_diag_scores = scores[source_landmark_idxs, target_landmark_idxs]
    
    # Obtain landmark-landmark off diag pairs
    L_L_off_diag_mask = np.logical_and(source_landmark_target_all_mask, source_all_target_landmark_mask)
    L_L_off_diag_mask &= ~landmark_landmark_mask
    L_L_off_diag_scores = scores[L_L_off_diag_mask]
    
    # Landmark - Non-landmark pairs
    L_non_L_mask = source_landmark_target_all_mask ^ source_all_target_landmark_mask
    L_non_L_scores = scores[L_non_L_mask]
    
    # Hom - Hom pairs
    H_H_mask = np.zeros_like(scores, dtype=bool)
    H_H_mask[source_homolog_idxs, target_homolog_idxs] = True
    H_H_mask[source_landmark_idxs, target_landmark_idxs] = False
    H_H_scores = scores[H_H_mask]
    
    # Obtain other scores
    other_mask = np.ones_like(scores, dtype=np.bool)
    other_mask &= ~(source_landmark_target_all_mask | source_all_target_landmark_mask | H_H_mask )
    other_scores = scores[other_mask] 

    return L_L_diag_scores,\
           L_L_off_diag_scores,\
           L_non_L_scores,\
           H_H_scores,\
            other_scores

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

def non_landmark_idxs(n, landmark_idxs):
    return [i for i in range(n) if i not in set(landmark_idxs)]

###############################################################################
# HANDL embedding
###############################################################################

def embed_matrices(source_C, target_D, landmark_idxs):
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

def embed_networks(source_G, target_G, homologs, n_landmarks, 
                 src_lam=0.05, tgt_lam=0.05, return_idxs=False):
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
    source_D = regularized_laplacian(source_G, source_nodes, src_lam)
    t_end = time.time()
    source_laplacian_time = t_end - t_start

    t_start = time.time()
    source_C = rkhs_factor(source_D)
    t_end = time.time()
    source_rkhs_time = t_end - t_start 

    t_start = time.time()
    target_D = regularized_laplacian(target_G, target_nodes, tgt_lam)
    t_end = time.time()
    target_laplacian_time = t_end - t_start

    source_node2index = util.node_to_index(source_G)
    target_node2index = util.node_to_index(target_G)

    landmark_homs = homologs[:n_landmarks]
    landmark_idxs = [(source_node2index[s_n], target_node2index[t_n]) for 
                     s_n, t_n in landmark_homs]

    t_start = time.time()
    source_C, target_C_hat = embed_matrices(source_C, target_D, landmark_idxs)
    t_end = time.time()
    handl_time = t_end - t_start

    runtimes = dict(source_regularized_laplacian_runtime=source_laplacian_time,
                    source_rkhs_factorization_runtime=source_rkhs_time,
                    target_regularized_laplacian_runtime=target_laplacian_time,
                    embed_runtime=handl_time)

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

def save_embeddings(X, nodes, landmarks, fp, weight=None):
    ''' Saves embedding data to given file '''
    joblib.dump(dict(X=X, nodes=nodes, landmarks=landmarks, weight=weight), fp)

#def homologs_in_graphs(G1, G2, homologs):
#    '''
#    Computes list of homologs from a pair of graphs and a given list of
#    homolog candidates
#    '''
#    G1_nodes = set(G1.nodes())
#    G2_nodes = set(G2.nodes())
#    valid_homologs = [(h1, h2) for h1, h2 in homologs \
#                        if h1 in G1_nodes and h2 in G2_nodes]
#    return valid_homologs
