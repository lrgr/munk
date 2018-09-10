import argparse
import os, sys

# Supress matplotlib display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from sklearn import neighbors
from sklearn.externals import joblib

import munk

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--munk_scores', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-xm', '--xmax', type=float, required=False, default=0.6)
    parser.add_argument('-lw', '--line_width', type=float, required=False, default=2)
    parser.add_argument('-fs', '--font_size', type=float, required=False, default=12)
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
    plt.tight_layout()
    plt.savefig(output)

def main(args):
    munk_data = joblib.load(args.munk_scores)
    homologs = munk_data['homologs']
    landmarks = munk_data['landmarks']
    A_nodes = munk_data['A_nodes']
    B_nodes = munk_data['B_nodes']
    sim_scores = munk_data['X']

    A_n2i = dict((n, i) for i, n in enumerate(A_nodes))
    B_n2i = dict((n, i) for i, n in enumerate(B_nodes))

    landmark_idxs = [(A_n2i[a], B_n2i[b]) for a, b in landmarks]
    homolog_idxs = [(A_n2i[a], B_n2i[b]) for a, b in homologs]

    dissims = 1. / sim_scores
    dissims /= np.mean(dissims)

    # Separate homolog, and other pairs from landmark pairs.
    _, _, _, hom_dissims, other_dissims = \
        munk.separate_scores(dissims, landmark_idxs, homolog_idxs)
    plots = [(hom_dissims, 'Homolog pairs'), (other_dissims, 'Other pairs')]
    plot_and_save(plots, 'Dissimilarity scores', args.output,
                 xmax=args.xmax,
                 font_size=args.font_size,
                 line_width=args.line_width)

if __name__ == '__main__':
    main(parse_args())
