#!/usr/bin/env python

# Load required modules
from collections import defaultdict

# Constants
SL_NAME = "SL"
NON_SL_NAME = "Non-SL"
UNKNOWN_NAME = "Inconclusive"
nameToIndex = dict(zip([NON_SL_NAME, SL_NAME, UNKNOWN_NAME], range(3)))

# Classify a gene as SL or not
def classify_gi(score, sl_threshold, uncertainty_threshold, pval, pval_threshold):
    score = float(score)
    if float(score) < sl_threshold and pval < pval_threshold:
        return SL_NAME
    # i.e. uncertainty_threshold < score < sl_threshold
    elif float(score) < uncertainty_threshold and pval < pval_threshold:
        return UNKNOWN_NAME
    else:
        return NON_SL_NAME

# Create a GI dataframe, categorizing pairs as SL or not.
# By default, we only require a score threshold. We optionally
# will use a p-value threshold as well, but p-values aren't given
# we just allow all pairs below the sl threshold.
def construct_gi_dataframe(pairToScore, sl_threshold, uncertainty_threshold, pairToPval=defaultdict(int), pval_threshold=1):
    import pandas as pd
    items = []
    for pair, score in pairToScore.iteritems():
        geneA, geneB = sorted(pair)
        items.append({
            "Gene A": geneA,
            "Gene B": geneB,
            "Score": float(score),
            "Category": classify_gi(score, sl_threshold, uncertainty_threshold,
                                    pairToPval[pair], pval_threshold)
        })

    df = pd.DataFrame(items)[['Gene A', 'Gene B', 'Score', 'Category']]
    return df.sort_values(['Gene A', 'Gene B'])

# Output a GI dataframe to file
def output_gi_dataframe(df, output_prefix, output_format):
    import numpy as np
    if output_format == 'tsv':
        df.to_csv(output_prefix + '.tsv', index=False, sep='\t')
    elif output_format == 'npy':
        # Create Numpy-style data
        genes = sorted(set(df['Gene A']) | set(df['Gene B']))
        n_genes = len(genes)
        geneToIndex = dict(zip(genes, range(n_genes)))

        A = [ [geneToIndex[r['Gene A']], geneToIndex[r['Gene B']], r['Score'], nameToIndex[r['Category']]]
              for i, r in df.iterrows() ]

        # Save the Numpy array and ordered list of nodes
        np.save(output_prefix + '.npy', A)
        with open(output_prefix + '-gene-index.txt', 'w') as OUT:
            OUT.write('\n'.join(genes))
    else:
        raise NotImplementedError('Output format "%s" not implemented.' % output_format)
