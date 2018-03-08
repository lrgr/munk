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
parser.add_argument('-s2t', '--source_to_target', type=str, default=None, required=False)
parser.add_argument('-nozero', '--nozero', type=str, default=None, required=False)
parser.add_argument('-renorm', '--renorm', type=str, default=None, required=False)
parser.add_argument('-v', '--verbosity', type=int, default=logging.INFO, required=False)

args = parser.parse_args(sys.argv[1:])

# Set up logger
logger = get_logger(args.verbosity)

################################################################################
# LOAD INPUT DATA
################################################################################
# Load features
logger.info('[Loading features and genetic interactions]')
src_data = joblib.load(args.feature_files[0])
tgt_data = joblib.load(args.feature_files[1])

X_S, y_S = np.array(src_data.get('X')), src_data.get('y')
S_pairs, S_name = src_data.get('pairs'), args.names[0]
S_args = src_data.get('args')
X_T, y_T = np.array(tgt_data.get('X')), tgt_data.get('y')
T_pairs, T_name = tgt_data.get('pairs'), args.names[1]
T_args = tgt_data.get('args')

from sklearn import preprocessing
import scipy as spy

logger.info('[Starting processing]')

if args.nozero:
    # normalize each feature to unit norm
    # note that assumptions of coral don't apply in this case
    NC_S = preprocessing.normalize(X_S, axis=0)
    NC_T = preprocessing.normalize(X_T, axis=0)
else:
    # normalize each feature to zero mean and unit norm
    NC_S = preprocessing.scale(X_S, axis=0)
    NC_T = preprocessing.scale(X_T, axis=0)

# compute covariances
CC_S = np.cov(NC_S.T)
CC_T = np.cov(NC_T.T)

# regularize (one could add a parameter here)
CS = CC_S + np.eye(*CC_S.shape)
CT = CC_T + np.eye(*CC_T.shape)

if args.source_to_target:
    # this is the recommended approach in the coral paper

    logger.info('[performing source to target fitting]')
    CS_sqrt = np.real_if_close(spy.linalg.sqrtm(CS))
    DS = NC_S.dot(np.linalg.pinv(CS_sqrt))
    NC_S = DS.dot(np.real_if_close(spy.linalg.sqrtm(CT)))

else:
    # contrary to the coral paper we are 'fitting' the target distribution
    # into the source's distributional shape

    logger.info('[performing target to source fitting]')
    CT_sqrt = np.real_if_close(spy.linalg.sqrtm(CT))
    DT = NC_T.dot(np.linalg.pinv(CT_sqrt))
    NC_T = DT.dot(np.real_if_close(spy.linalg.sqrtm(CS)))

if args.renorm:
    NC_S = preprocessing.normalize(NC_S, axis=0)
    NC_T = preprocessing.normalize(NC_T, axis=0)

# Output to files
# Rewriting existing files for now so I don't have to change the snakefile
output = dict(X=NC_S, y=y_S, pairs=S_pairs, params=S_args)
joblib.dump(output, args.feature_files[0])
output = dict(X=NC_T, y=y_T, pairs=T_pairs, params=T_args)
joblib.dump(output, args.feature_files[1])
logger.info('[Done]')

    
