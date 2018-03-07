#!/usr/bin/env python

################################################################################
# SETUP AND PARSING COMMAND-LINE ARGUMENTS
################################################################################
# Load required modules
import sys, os, argparse, numpy as np, logging, pandas as pd
from sklearn.externals import joblib
from copy import deepcopy

# Load our modules
from constants import *
this_dir = os.path.dirname(__file__)
sys.path.append(os.path.join(this_dir, '../../src'))
from i_o import get_logger

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-ff', '--feature_files', type=str, required=True, nargs=2)
parser.add_argument('-n', '--names', type=str, required=True, nargs=2)
parser.add_argument('-o', '--output_file', type=str, required=True)
parser.add_argument('-v', '--verbosity', type=int, default=logging.INFO, required=False)
parser.add_argument('-rs', '--random_seed', type=int, default=1764, required=False)
parser.add_argument('-nf', '--n_folds', type=int, default=4, required=False)
parser.add_argument('-nj', '--n_jobs', type=int, default=1, required=False)
parser.add_argument('-ho', '--hold-out', type=str, choices=[GENE_PAIRS, GENES],
                    default=GENE_PAIRS, required=False)

# Add classifier choices with subparsers
subparser = parser.add_subparsers(dest='classifier', help='Classifier')
rf_parser = subparser.add_parser('rf')
rf_parser.add_argument('-md', '--max_depth', type=int, default=None, required=False)
rf_parser.add_argument('-nt', '--n_trees', type=int, default=100, required=False)

svm_parser = subparser.add_parser('svm')
svm_parser.add_argument('-sc', '--svm_C', type=float, default=1.0, required=False)
svm_parser.add_argument('-st', '--svm_tolerance', type=float, default=1e-3, required=False)

args = parser.parse_args(sys.argv[1:])

# Set up logger
logger = get_logger(args.verbosity)

################################################################################
# LOAD INPUT DATA
################################################################################
# Load features
logger.info('[Loading features and genetic interactions]')
a_data = joblib.load(args.feature_files[0])
b_data = joblib.load(args.feature_files[1])

X_A, y_A = np.array(a_data.get('X')), a_data.get('y')
A_pairs, A_name = a_data.get('pairs'), args.names[0]
X_B, y_B = np.array(b_data.get('X')), b_data.get('y')
B_pairs, B_name = b_data.get('pairs'), args.names[1]

# Log info about the data
logger.info('- Species A (%s): %s samples x %s features (%s SLs)' % (A_name, X_A.shape[0], X_A.shape[1], int(y_A.sum())))
logger.info('- Species B (%s): %s samples x %s features (%s SLs)' % (B_name, X_B.shape[0], X_B.shape[1], int(y_B.sum())))

################################################################################
# TRAIN WITHIN ONE SPECIES AND PREDICT IN THE OTHER
################################################################################
# Load the required modules
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.model_selection import KFold
from sklearn.metrics import average_precision_score, roc_auc_score, precision_recall_curve
from itertools import permutations

# Simple function for outputting results
def format_result(results, skip_keys={"Name"}):
    return '; '.join('{0}={1:g}'.format(k, v) for k, v in results.items() if k not in skip_keys)

# Max F1 score
def f1(p, r):
    if p == 0 and r == 0:
        return 0.
    else:
        return 2.0 * p * r / (p + r)
    
def max_f1_score(_y_true, _y_hat):
    ps, rs, _ = precision_recall_curve(_y_true, _y_hat)
    f1s = np.asarray([ f1(p, r) for p, r in zip(ps, rs)])
    return max(f1s)

# Helper functions to be used by Parallel to train/predict
def data_producer(random_state):
    kf = KFold(n_splits=args.n_folds, random_state=random_state, shuffle=True)
    if args.hold_out == GENE_PAIRS:
        for i, (train_index, test_index) in enumerate(kf.split(X_A)):
            # Split up the source data
            X_train, X_test = X_A[train_index], X_A[test_index]
            y_train, y_test = y_A[train_index], y_A[test_index]

            # Train the classifier within the source and predict within and across species
            random_state += 1

            yield i+1, X_train, y_train, X_test, y_test, random_state
    elif args.hold_out == GENES:
        # Get a list of genes
        A_genes = sorted(set( g for p in A_pairs for g in p ))

        # Split on genes
        A_pair_indices = set(range(len(A_pairs)))
        for i, (train_gene_index, test_gene_index) in enumerate(kf.split(A_genes)):
            # Get the train/test indices
            train_gene_set = set( A_genes[i] for i in train_gene_index )
            train_index = [ j for j, p in enumerate(A_pairs) if all( g in train_gene_set for g in p ) ]
            test_index = sorted(A_pair_indices - set(train_index))

            # Split up the source data
            X_train, X_test = X_A[train_index], X_A[test_index]
            y_train, y_test = y_A[train_index], y_A[test_index]

            # Increment random state
            random_state += 1
            yield i+1, X_train, y_train, X_test, y_test, random_state

def train_and_predict(fold, _X_train, _y_train, _X_test, _y_test, _random_state):
    # Random forest
    if args.classifier == 'rf':
        # Train the random forest
        clf = RandomForestClassifier(n_estimators=args.n_trees, max_depth=args.max_depth,
                                     random_state=_random_state,
                                     n_jobs=args.n_jobs)
        clf.fit(_X_train, _y_train)

        # Make out-of-sample predictions
        y_test_hat = clf.predict_proba(_X_test)[:, 1]
        y_B_hat = clf.predict_proba(X_B)[:, 1]
        
    # SVM
    elif args.classifier == 'svm':
        # Train the Linear SVM
        clf = LinearSVC(C=args.svm_C, tol=args.svm_tolerance)
        clf.fit(_X_train, _y_train)

        # Make out of sample predictions. Decision function outputs
        # a single list of values for binary class data.
        y_test_hat = clf.decision_function(_X_test)
        y_B_hat = clf.decision_function(X_B)
        
    else:
        raise NotImplementedError('Classifier "%s" not implemented.' % args.classifier)

        
    # Report the results
    result = {
        "Source": {
            "Name": A_name,
            "AUPRC": average_precision_score(_y_test, y_test_hat),
            "AUROC": roc_auc_score(_y_test, y_test_hat),
            "F1": max_f1_score(_y_test, y_test_hat)
        },
        "Target": {
            "Name": B_name,
            "AUPRC": average_precision_score(y_B, y_B_hat),
            "AUROC": roc_auc_score(y_B, y_B_hat),
            "F1": max_f1_score(y_B, y_B_hat)
        },
        "Fold": fold,
        "Train size": len(_y_train),
        "Test size": len(_y_test),
    }
    logger.info('- Fold: {}'.format(fold))
    logger.info('\t- {0}->{0}: {1}'.format(A_name, format_result(result['Source'])))
    logger.info('\t- {0}->{1}: {2}'.format(A_name, B_name, format_result(result['Target'])))
    
    return result

# Train on A, predict on held-out A and B, executing in parallel
logger.info('[Training and evaluating models]')
from sklearn.externals.joblib import Parallel, delayed
r = Parallel(n_jobs=min(args.n_folds, args.n_jobs), verbose=0)( delayed(train_and_predict)(*d) for d in data_producer(args.random_seed) )
results = r

# Add an "average" row across folds, and report the current results
average_result = deepcopy(results[-1])
average_result['Fold'] = 'Average'
for s_name in ['Source', 'Target']:
    for measure in ['AUPRC', 'AUROC', 'F1']:
        average_result[s_name][measure] = np.mean([ r[s_name][measure] for r in results ])
        
results.append(average_result)

# Flatten results
flat_results = []
for r in results:
    for ty in ['Source', 'Target']:
        flat_results.append({
            "Train": r['Source']['Name'],
            "Test": r[ty]['Name'],
            "AUROC": r[ty]['AUROC'],
            "AUPRC": r[ty]['AUPRC'],
            "F1": r[ty]['F1'],
            "Fold": r['Fold'],
            "Train size": r['Train size'],
            "Test size": r['Test size']
        })
    
# Output results to file
df = pd.DataFrame(flat_results)[['Train', 'Test', 'Fold', 'F1', 'AUROC', 'AUPRC', 'Train size', 'Test size']]
df = df.sort_values(['Train', 'Test', 'Fold'])
df.to_csv(args.output_file, sep='\t', index=False)
