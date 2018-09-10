import argparse
import json
import random
import time

import networkx as nx
import numpy as np
from sklearn.externals import joblib

import munk
import munk.util as util

###############################################################################
# Arguments
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
    parser.add_argument('-rs', '--random_seed', type=int, required=False, default=28791)
    parser.add_argument('--src_lam', type=float, required=False, default=0.05)
    parser.add_argument('--tgt_lam', type=float, required=False, default=0.05)
    return parser.parse_args()

###############################################################################
# MAIN
###############################################################################

def main(args):
    '''
    Compute MUNK embeddings given target and species edgelists and list of
    homologs.
    '''
    log = munk.io.get_logger()
    random.seed(args.random_seed)

    log.info('Loading homologs list from %s', args.homolog_list)
    raw_homologs = util.read_homolog_list(args.homolog_list)

    log.info('Loading source edgelist from %s', args.source_edgelist)
    source_G = nx.read_edgelist(args.source_edgelist, encoding='ascii')
    log.info('Loading target edgelist from %s', args.target_edgelist)
    target_G = nx.read_edgelist(args.target_edgelist, encoding='ascii')

    log.info('Computing MUNK embeddings with %d landmarks', args.n_landmarks)
    n_landmarks = args.n_landmarks

    source_G = util.simple_two_core(source_G)
    target_G = util.simple_two_core(target_G)
    homologs = util.homologs_in_graphs(source_G, target_G, raw_homologs)

    t_start = time.time()
    source_data, target_data, landmarks, runtimes = \
        munk.embed_networks(source_G, target_G, homologs, n_landmarks,
                             src_lam=args.src_lam, tgt_lam=args.tgt_lam)
    t_end = time.time()
    total_time = t_end - t_start
    source_X, source_nodes = source_data
    target_X, target_nodes = target_data
    source_landmarks = [ l[0] for l in landmarks ]
    target_landmarks = [ l[1] for l in landmarks ]

    log.info('Saving source species embeddings to %s', args.source_output_file)
    munk.save_embeddings(source_X, source_nodes, source_landmarks, args.source_output_file)
    log.info('Saving target species embeddings to %s', args.target_output_file)
    munk.save_embeddings(target_X, target_nodes, target_landmarks, args.target_output_file)

    log.info('Source data shape {}'.format(source_X.shape))
    log.info('Target data shape {}'.format(target_X.shape))

    log.info('Saving MUNK similarity scores to %s', args.sim_scores_output_file)

    munk_scores = np.dot(source_X, target_X.T)
    joblib.dump(dict(X=munk_scores, A_nodes=source_nodes,
                     B_nodes=target_nodes, landmarks=landmarks,
                     homologs=homologs),
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

    log.info('MUNK embedding complete!')

if __name__ == '__main__':
    main(parse_args())
