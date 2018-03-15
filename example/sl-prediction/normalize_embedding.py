#!/usr/bin/env python

# Load required modules
import sys, os, logging, numpy as np
this_dir = os.path.dirname(__file__)
sys.path.append(os.path.join(this_dir, '../../src'))

# CORAL algorithm (Sun, Feng, & Saenko, arXiv:1612.01939v1)
def coral(D_S, D_T, source_to_target=True, lamb=0.00001):
    """
    Arguments:
    - D_S: Source data (numpy.array)
    - D_T: Target data (numpy.array)
    - source_to_target: Transform source into target space (bool)
    Returns:
    - D_S^*: Adjusted source data
    """
    # Set up
    from sklearn.preprocessing import scale
    from scipy.linalg import sqrtm
    from numpy.linalg import pinv, inv
    C_S = np.cov(D_S.T) + lamb*np.eye(D_S.shape[1])
    C_T = np.cov(D_T.T) + lamb*np.eye(D_T.shape[1])
    C_S_centered = scale(D_S, with_std = False)
    C_T_centered = scale(D_T, with_std = False)
    
    # Perform CORAL (Algorithm 1)
    if source_to_target:
        # this is the recommended approach in the coral paper
        C_S_sqrt = np.real_if_close(sqrtm(C_S))
        D_S = C_S_centered @ pinv(C_S_sqrt)
        return D_S @ np.real_if_close(sqrtm(C_T))
    else:
        # contrary to the coral paper we are 'fitting' the target distribution
        # into the source's distributional shape
        C_T_sqrt = np.real_if_close(sqrtm(C_T))
        D_T = C_T_centered @ inv(C_T_sqrt)
        return D_T @ np.real_if_close(sp.linalg.sqrtm(C_S))

if __name__ == '__main__':
    # Parse command line arguments
    import argparse
    from constants import *
    parser = argparse.ArgumentParser()
    parser.add_argument('-sf', '--source_embedding_file', type=str, required=True)
    parser.add_argument('-tf', '--target_embedding_file', type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-n', '--normalization', type=str, required=True,
                        choices=[CORAL_NORM, NONE_NORM])
    parser.add_argument('-v', '--verbosity', type=int, default=logging.INFO, required=False)
    args = parser.parse_args(sys.argv[1:])

    from i_o import get_logger
    logger = get_logger(args.verbosity)

    # Load the embeddings
    from sklearn.externals import joblib
    logger.info('[Loading embeddings]')
    obj_A = joblib.load(args.source_embedding_file)
    X_A = obj_A.get('X')
    nodes_A = obj_A.get('nodes')
    landmarks_A = obj_A.get('landmarks')
    weight_A = obj_A.get('weight')
    logger.info('- Source embedding: %s nodes' % len(nodes_A))

    obj_B = joblib.load(args.target_embedding_file)
    X_B = obj_B.get('X')
    nodes_B = obj_B.get('nodes')
    logger.info('- Target embedding: %s nodes' % len(nodes_A))

    assert( X_A.shape[1] == X_B.shape[1] )

    # Apply normalization
    if args.normalization == CORAL_NORM:
        logger.info('[Aligning feature spaces with CORAL]')
        X_A = coral(X_A, X_B)
    elif args.normalization == NONE_NORM:
        logger.info('[Peforming no normalization]')
    else:
        raise NotImplementedError('Normalization mode "%s" not implemented.' % args.normalization)

    # Output aligned source embedding to file
    from util import save_embeddings
    save_embeddings(X_A, nodes_A, landmarks_A, args.output_file, weight_A)
    
