from __future__ import division
import numpy as np
import argparse
import scipy
import scipy.stats
import scipy.linalg
import sklearn
import sklearn.neighbors
import time

from sklearn.externals import joblib
import logging_utils

import matplotlib
import os
import warnings
with warnings.catch_warnings():
    # ignore the warning if matplotlib has already been loaded
    warnings.simplefilter('ignore')
    # so that access to an X server is not needed -
    # can run without a terminal
    # must be called before pyplot is imported
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-hf', '--homologs_file', type=str, required=True)
    parser.add_argument('-i', '--input_handl_file', type=str, required=True)
    parser.add_argument('-pl', '--plot_lim', type=float, default = 0.01)
    parser.add_argument('-o', '--output_file', type=str, default='homolog-density')
    args = parser.parse_args()

    data = joblib.load(args.input_handl_file)
    target_handl_C =              data['target_handl_C']
    source_C =                    data['source_C']
    D =                           data['D']
    target_nodes =                data['target_nodes']
    source_nodes =                data['source_nodes']
    target_landmarks =            data['target_landmarks']
    target_landmark_indices =     data['target_landmark_indices']
    target_non_landmark_indices = data['target_non_landmark_indices']
    source_landmarks =            data['source_landmarks']
    source_landmark_indices =     data['source_landmark_indices']
    source_non_landmark_indices = data['source_non_landmark_indices']

    # for this plot we want to work with dissimilarities
    D = 1./D

    # normalize
    D = D / np.mean(D)

    # read in homologs
    # assume the order is source node - target node
    source_homs = []
    target_homs = []
    with open(args.homologs_file, 'r') as hom_f:
        for line in hom_f:
            source_hom, target_hom = line.split()
            if (source_hom in source_nodes) and (target_hom in target_nodes):
                source_homs.append(source_hom)
                target_homs.append(target_hom)

    # we don't want to use known landmark-homologs
    # so we remove those entries from the HANDL score matrix
    source_valid = [r for r in range(D.shape[0])
                        if r not in source_landmark_indices]
    target_valid = [r for r in range(D.shape[1])
                        if r not in target_landmark_indices]
    D = D[source_valid][:,target_valid]

    # update the names of rows and columns
    source_nodes = [n for i, n in enumerate(source_nodes) if i in source_valid]
    target_nodes = [n for i, n in enumerate(target_nodes) if i in target_valid]

    # update these indices as well
    source_non_landmark_indices = [source_nodes.index(node) for node in source_homs
                                       if node not in source_landmarks]
    target_non_landmark_indices = [target_nodes.index(node) for node in target_homs
                                       if node not in target_landmarks]
    
    # get HANDL scores for Non-landmark homologs
    NLH = D[source_non_landmark_indices,target_non_landmark_indices]

    # create kernel density estimator
    kde = sklearn.neighbors.KernelDensity(kernel='gaussian',
                                              bandwidth = args.plot_lim / 20.)

    # need to add another dimension as required by sklearn
    # arrays passed to kde must be 2-dimensional
    X_plot = np.reshape(np.linspace(0, args.plot_lim, 500), (-1, 1))

    # Density of non-landmark homolog pairs
    kde.fit(np.reshape(np.ravel(NLH), (-1, 1)))
    NLHdens = kde.score_samples(X_plot)

    # Density of all pairs
    Dsamp = np.random.choice(np.ravel(D),100000)
    kde.fit(np.reshape(Dsamp, (-1, 1)))
    Ddens = kde.score_samples(X_plot)

    # create and save plot
    font_size = 20
    line_width = 3
    plt.figure()
    plt.plot(X_plot[:,0], np.exp(NLHdens), lw = line_width,
                 label = 'Homolog Pairs')
    plt.plot(X_plot[:,0], np.exp(Ddens), '-.',
                 lw = line_width, ls = (0,(5,2)), label = 'All Pairs')
    plt.ylabel('Density', size = font_size)
    plt.xlabel('HANDL dissimilarity', size = font_size)
    plt.legend(loc='best', fontsize = font_size)
    plt.savefig('{}.pdf'.format(args.output_file))
    plt.close()
