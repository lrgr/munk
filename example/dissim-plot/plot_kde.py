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
    parser.add_argument('-i', '--HANDL_scores', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('--xmax', type=float, required=False, default=0.6)
    return parser.parse_args()
    

def plot_and_save(scores_and_labels, xlabel, output,
                  xmin=0.0, 
                  xmax=0.6, 
                  smoothness=20., 
                  font_size=12, 
                  line_width=2):
    
    # create and save plot
    plt.figure()
    
    # create kernel density estimator
    kde = neighbors.KernelDensity(kernel='gaussian', bandwidth = xmax / smoothness)
    # need to add another dimension as required by sklearn
    # arrays passed to kde must be 2-dimensional
    X_plot = np.reshape(np.linspace(xmin, xmax, 500), (-1, 1))
    styles = ['-', '--', '-.', ':']
    for i, (xs, label) in enumerate(scores_and_labels):
        scores = np.ravel(xs) if len(xs) < 1e5 else np.random.choice(np.ravel(xs), int(1e5))
        kde.fit(np.reshape(scores, (-1, 1)))
        densities = kde.score_samples(X_plot)
        plt.plot(X_plot[:,0], np.exp(densities), lw = line_width,
                 label = label, ls=styles[i % len(styles)])
    plt.ylabel('Density', size = font_size)
    plt.xlabel(xlabel, size = font_size)
    plt.legend(loc='best', fontsize = font_size)
    plt.savefig(output)

def main(args):
    HANDL_data = joblib.load(args.HANDL_scores)
    homologs = HANDL_data['homologs']
    landmarks = HANDL_data['landmarks']
    A_nodes = HANDL_data['A_nodes']
    B_nodes = HANDL_data['B_nodes']
    sim_scores = HANDL_data['X']

    A_n2i = dict((n, i) for i, n in enumerate(A_nodes))
    B_n2i = dict((n, i) for i, n in enumerate(B_nodes))

    landmark_idxs = [(A_n2i[a], B_n2i[b]) for a, b in landmarks]
    homolog_idxs = [(A_n2i[a], B_n2i[b]) for a, b in homologs]
    
    dissims = 1. / sim_scores
    dissims /= np.mean(dissims)
    _, _, _, hom_dissims, other_dissims = \
        handl.separate_scores(dissims, landmark_idxs, homolog_idxs)
    plots = [(hom_dissims, 'Homologs'), (other_dissims, 'Other')]
    plot_and_save(plots, 'Dissimilarity scores', args.output)

    
    # Load 
    pass

if __name__ == '__main__':
    main(parse_args())
