import argparse
import os, sys

import matplotlib

# Supress matplotlib display
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import numpy as np
from sklearn import neighbors
from sklearn.externals import joblib

sys.path.append("../../src")
import handl
import i_o

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hs', '--HANDL_scores', required=True)
    parser.add_argument('-rs', '--resnik_scores', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-lw', '--line_width', type=float, required=False, default=2)
    parser.add_argument('-fs', '--font_size', type=float, required=False, default=12)
    return parser.parse_args()

def main(args):
    HANDL_data = joblib.load(args.HANDL_scores)
    homologs = HANDL_data['homologs']
    landmarks = HANDL_data['landmarks']
    A_nodes = HANDL_data['A_nodes']
    B_nodes = HANDL_data['B_nodes']
    HANDL_scores_raw = HANDL_data['X']

    A_n2i = dict((n, i) for i, n in enumerate(A_nodes))
    B_n2i = dict((n, i) for i, n in enumerate(B_nodes))

    # Get Resnik scores 
    resnik_data = np.load(args.resnik_scores)
    resnik_A_n2i = resnik_data['leftIndex'][()]
    resnik_B_n2i = resnik_data['rightIndex'][()]
    resnik_scores_raw = resnik_data['Rscore']

    print('# A genes with scores:', len(set(A_nodes) & set(resnik_A_n2i.keys())))
    print('# B genes with scores:', len(set(B_nodes) & set(resnik_B_n2i.keys())))

    A_landmarks, B_landmarks = [set(ls) for ls in zip(*landmarks)]

    A_scored_nodes = sorted((set(A_nodes) & set(resnik_A_n2i.keys())) - A_landmarks)
    B_scored_nodes = sorted((set(B_nodes) & set(resnik_B_n2i.keys())) - B_landmarks)

    print('# scored pairs:', len(A_scored_nodes) * len(B_scored_nodes))
    
    HANDL_scores = []
    resnik_scores = []
    r_A_idxs = [resnik_A_n2i[node] for node in A_scored_nodes]
    r_B_idxs = [resnik_B_n2i[node] for node in B_scored_nodes]

    h_A_idxs = [A_n2i[node] for node in A_scored_nodes]
    h_B_idxs = [B_n2i[node] for node in B_scored_nodes]

    for r_A, h_A in zip(r_A_idxs, h_A_idxs):
        for r_B, h_B in zip(r_B_idxs, h_B_idxs):
            HANDL_scores.append(HANDL_scores_raw[h_A, h_B])
            resnik_scores.append(resnik_scores_raw[r_A, r_B])

    HANDL_scores = np.asarray(HANDL_scores)
    resnik_scores = np.asarray(resnik_scores)
    sort_idxs = np.argsort(HANDL_scores)[::-1]
    n_scores = len(sort_idxs)
    rand_idxs =  np.random.permutation(n_scores)
    resnik_scores = np.take(resnik_scores, sort_idxs)
    rand_resnik_scores = np.take(resnik_scores, rand_idxs)

    #smooth over bins
    binsize = 100000
    binned_scores = n_scores - (n_scores % binsize)
    resnik_scores = np.nanmean(resnik_scores[:binned_scores].reshape((-1,binsize)), axis=1)
    rand_resnik_scores = np.nanmean(rand_resnik_scores[:binned_scores].reshape((-1,binsize)), axis=1)
    n_bins = len(resnik_scores)
    print('# bins', n_bins)
    plt.figure()
    plt.plot(np.arange(n_bins), rand_resnik_scores, 
             label='Ranked randomly', lw = args.line_width)
    plt.plot(np.arange(n_bins), resnik_scores,
             label='Ranked by HANDL similarity', lw=args.line_width)
    plt.xlabel('Ranked pairs',size=args.font_size)
    plt.ylabel('Resnik score', size=args.font_size)
    plt.legend(loc='best', fontsize=args.font_size)
    plt.savefig(args.output)
    plt.close()



if __name__ == '__main__':
    main(parse_args())
