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
parser.add_argument('-nif', '--n_inner_folds', type=int, default=3, required=False)
parser.add_argument('-nj', '--n_jobs', type=int, default=1, required=False)
parser.add_argument('-ho', '--hold-out', type=str, choices=[GENE_PAIRS, GENES],
                    default=GENE_PAIRS, required=False)
#TODO 
#parser.add_argument('-po', '--predictions_output_file', type=str, required=True)

# Add classifier choices with subparsers
subparser = parser.add_subparsers(dest='classifier', help='Classifier')
rf_parser = subparser.add_parser('rf')
rf_parser.add_argument('-md', '--max_depth', type=int, default=None, required=False)
rf_parser.add_argument('-nt', '--n_trees', type=int, nargs='*', required=False,
                       default=[10, 100, 250, 500])

svm_parser = subparser.add_parser('svm')
svm_parser.add_argument('-sc', '--svm_Cs', type=float, required=False, nargs='*',
                        default=[0.01, 0.1, 1, 10, 100, 1000, 10000])
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
A_pairs, A_name = np.asarray(a_data.get('pairs')), args.names[0]
X_B, y_B = np.array(b_data.get('X')), b_data.get('y')
B_pairs, B_name = np.asarray(b_data.get('pairs')), args.names[1]

# Log info about the data
logger.info('- Species A (%s): %s samples x %s features (%s SLs)' % (A_name, X_A.shape[0], X_A.shape[1], int(y_A.sum())))
logger.info('- Species B (%s): %s samples x %s features (%s SLs)' % (B_name, X_B.shape[0], X_B.shape[1], int(y_B.sum())))

################################################################################
# TRAIN WITHIN ONE SPECIES AND PREDICT IN THE OTHER
################################################################################
# Load the required modules
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.model_selection import KFold, GridSearchCV
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
        for i, ((A_train_index, A_test_index), (B_train_index, B_test_index)) in enumerate(zip(kf.split(X_A), kf.split(X_B))):

            # Split up the source datia
            A_pairs_train, A_pairs_test = A_pairs[A_train_index], A_pairs[A_test_index]
            X_A_train, X_A_test = X_A[A_train_index], X_A[A_test_index]
            y_A_train, y_A_test = y_A[A_train_index], y_A[A_test_index]
            A_data = (X_A_train, X_A_test, y_A_train, y_A_test, A_pairs_train, A_pairs_test)

            # Split up the train data
            B_pairs_train, B_pairs_test = B_pairs[B_train_index], B_pairs[B_test_index]
            X_B_train, X_B_test = X_B[B_train_index], X_B[B_test_index]
            y_B_train, y_B_test = y_B[B_train_index], y_B[B_test_index]
            B_data = (X_B_train, X_B_test, y_B_train, y_B_test, B_pairs_train, B_pairs_test)
            # Train the classifier within the source and predict within and across species
            random_state += 1

            yield i+1, A_data, B_data, random_state

    elif args.hold_out == GENES:
        # Get a list of genes
        A_genes = sorted(set( g for p in A_pairs for g in p ))
        B_genes = sorted(set( g for p in B_pairs for g in p ))

        # Split on genes
        A_pair_indices = set(range(len(A_pairs)))
        B_pair_indices = set(range(len(B_pairs)))
        for i, ((A_train_gene_index, A_test_gene_index), (B_train_gene_index, B_test_gene_index)) in enumerate(zip(kf.split(A_genes), kf.split(B_genes))):
            # Get the train/test indices
            A_train_gene_set = set( A_genes[i] for i in A_train_gene_index )
            A_train_index = [ j for j, p in enumerate(A_pairs) if all( g in A_train_gene_set for g in p ) ]
            A_test_index = sorted(A_pair_indices - set(A_train_index))

            B_train_gene_set = set( B_genes[i] for i in B_train_gene_index )
            B_train_index = [ j for j, p in enumerate(B_pairs) if all( g in B_train_gene_set for g in p ) ]
            B_test_index = sorted(B_pair_indices - set(B_train_index))

            # Split up the source datia
            A_pairs_train, A_pairs_test = A_pairs[A_train_index], A_pairs[A_test_index]
            X_A_train, X_A_test = X_A[A_train_index], X_A[A_test_index]
            y_A_train, y_A_test = y_A[A_train_index], y_A[A_test_index]
            A_data = (X_A_train, X_A_test, y_A_train, y_A_test, A_pairs_train, A_pairs_test)

            # Split up the train data
            B_pairs_train, B_pairs_test = B_pairs[B_train_index], B_pairs[B_test_index]
            X_B_train, X_B_test = X_B[B_train_index], X_B[B_test_index]
            y_B_train, y_B_test = y_B[B_train_index], y_B[B_test_index]
            B_data = (X_B_train, X_B_test, y_B_train, y_B_test, B_pairs_train, B_pairs_test)
            # Train the classifier within the source and predict within and across species
            random_state += 1

            yield i+1, A_data, B_data, random_state


def train_and_predict(fold, _A_data, _B_data, _random_state):
    inner_cv = KFold(n_splits=args.n_inner_folds, random_state=_random_state, shuffle=True)
    # unpack data...
    _X_A_train, _X_A_test, _y_A_train, _y_A_test, _A_pairs_train, _A_pairs_test = _A_data
    _X_B_train, _X_B_test, _y_B_train, _y_B_test, _B_pairs_train, _B_pairs_test = _B_data
    
    # concatenate
    _X_train = np.concatenate((_X_A_train, _X_B_train))
    _y_train = np.concatenate((_y_A_train, _y_B_train))

    # Random forest
    if args.classifier == 'rf':
        # Train the random forest
        rf = RandomForestClassifier(n_estimators=args.n_trees, max_depth=args.max_depth,
                                     random_state=_random_state,
                                     n_jobs=args.n_jobs)
        clf = GridSearchCV(rf, dict(n_estimators=args.n_trees), cv=inner_cv, 
                           n_jobs=1, pre_dispatch=1,
                           refit=True, scoring='average_precision')
        clf.fit(_X_train, _y_train)
        best_params = clf.best_params_

        # Make out-of-sample predictions
        y_A_test_hat = clf.predict_proba(_X_A_test)[:, 1]
        y_B_test_hat = clf.predict_proba(_X_B_test)[:, 1]
        
    # SVM
    elif args.classifier == 'svm':
        # Train the Linear SVM
        svc = LinearSVC(tol=args.svm_tolerance, random_state=_random_state+1)
        clf = GridSearchCV(svc, dict(C=args.svm_Cs), cv=inner_cv, 
                           n_jobs=args.n_jobs, pre_dispatch=args.n_jobs,
                           refit=True, scoring='average_precision')
        clf.fit(_X_train, _y_train)
        best_params = clf.best_params_

        # Make out of sample predictions. Decision function outputs
        # a single list of values for binary class data.
        y_A_test_hat = clf.decision_function(_X_A_test)
        y_B_test_hat = clf.decision_function(_X_B_test)
    else:
        raise NotImplementedError('Classifier "%s" not implemented.' % args.classifier)
        
    # Report the results
    result = {
        "Source": {
            "Name": A_name,
            "AUPRC": average_precision_score(_y_A_test, y_A_test_hat),
            "AUROC": roc_auc_score(_y_A_test, y_A_test_hat),
            "F1": max_f1_score(_y_A_test, y_A_test_hat),
            "Test size": len(_y_A_test),
        },
        "Target": {
            "Name": B_name,
            "AUPRC": average_precision_score(_y_B_test, y_B_test_hat),
            "AUROC": roc_auc_score(_y_B_test, y_B_test_hat),
            "F1": max_f1_score(_y_B_test, y_B_test_hat),
            "Test size": len(_y_B_test),
        },
        "Fold": fold,
        "Train size": len(_y_train),
        "Best params": str(best_params)
    }
    logger.info('- Fold: {} ({})'.format(fold, str(best_params)))
    logger.info('\t- {0}: {1}'.format(A_name, format_result(result['Source'])))
    logger.info('\t- {0}: {1}'.format(B_name, format_result(result['Target'])))
    
    return result, None #(clf, _pairs_test, _y_test, y_test_hat, B_pairs, y_B, y_B_hat)

# Train on A, predict on held-out A and B, executing in parallel
logger.info('[Training and evaluating models]')
from sklearn.externals.joblib import Parallel, delayed
#r = Parallel(n_jobs=min(args.n_folds, args.n_jobs), verbose=0)( delayed(train_and_predict)(*d) for d in data_producer(args.random_seed) )
r = [train_and_predict(*d) for d in data_producer(args.random_seed)]
results, clfs_and_preds = zip(*r)
results = list(results)

# Add an "average" row across folds, and report the current results
average_result = deepcopy(results[-1])
average_result['Fold'] = 'Average'
average_result['Best params'] = 'N/A'
for s_name in ['Source', 'Target']:
    for measure in ['AUPRC', 'AUROC', 'F1']:
        average_result[s_name][measure] = np.mean([ r[s_name][measure] for r in results ])
        
results.append(average_result)

logger.info('- Average')
logger.info('\t- {0}: {1}'.format(A_name, format_result(average_result['Source'])))
logger.info('\t- {0}: {1}'.format(B_name, format_result(average_result['Target'])))

# Flatten results
flat_results = []
for r in results:
    for ty in ['Source', 'Target']:
        flat_results.append({
            "Train": '{}+{}'.format(r['Source']['Name'], r['Target']['Name']),
            "Test": r[ty]['Name'],
            "AUROC": r[ty]['AUROC'],
            "AUPRC": r[ty]['AUPRC'],
            "F1": r[ty]['F1'],
            "Fold": r['Fold'],
            "Train size": r['Train size'],
            "Test size": r[ty]['Test size'],
            "Best params": r['Best params']
        })
    
# Output results to file
df = pd.DataFrame(flat_results)[['Train', 'Test', 'Fold', 'F1', 'AUROC', 'AUPRC', 'Train size', 'Test size', 'Best params']]
df = df.sort_values(['Train', 'Test', 'Fold'])
df.to_csv(args.output_file, sep='\t', index=False)

# TODO
#results_data = dict(summary=flat_results, data=[])
#
#for clf, pairs_test, y_test, y_test_hat, B_pairs, y_B, y_B_hat in clfs_and_preds:
#    results_data['data'].append(
#        dict(clf_type=args.classifier,
#             clf=clf,
#             pairs_test=pairs_test,
#             y_test=y_test,
#             y_test_hat=y_test_hat,
#             B_pairs=B_pairs,
#             y_B=y_B,
#             y_B_hat=y_B_hat))
## Save SVM, and predictions to disk...
#joblib.dump(results_data, args.predictions_output_file)
