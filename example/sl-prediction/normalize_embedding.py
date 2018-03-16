#!/usr/bin/env python

# Load required modules
import sys, os, logging, numpy as np
this_dir = os.path.dirname(__file__)
sys.path.append(os.path.join(this_dir, '../../src'))

# CORAL algorithm (Sun, Feng, & Saenko, arXiv:1612.01939v1)
def coral(D_S, D_T, lamb=0.00001):
    """
    Arguments:
    - D_S: Source data (numpy.array)
    - D_T: Target data (numpy.array)
    - lamb: Regularization parameter (float)
    Returns:
    - D_S^*: Adjusted, centered source data
    - D_T: Centered target data
    """
    # Set up
    from sklearn.preprocessing import scale
    from scipy.linalg import sqrtm
    from numpy.linalg import pinv, inv
    C_S = np.cov(D_S.T) + lamb*np.eye(D_S.shape[1])
    C_T = np.cov(D_T.T) + lamb*np.eye(D_T.shape[1])
    D_S_centered = scale(D_S, with_std = False)
    D_T_centered = scale(D_T, with_std = False)
    
    # Perform CORAL (Algorithm 1)
    C_S_sqrt = np.real_if_close(sqrtm(C_S))
    D_S = D_S_centered @ pinv(C_S_sqrt)
    D_S_star = D_S @ np.real_if_close(sqrtm(C_T))
    
    # Then center and return
    D_S_star_centered = scale(D_S_star, with_std=False)
    
    return D_S_star_centered, D_T_centered

if __name__ == '__main__':
    # Parse command line arguments
    import argparse
    from constants import *
    parser = argparse.ArgumentParser()
    parser.add_argument('-sf', '--source_embedding_file', type=str, required=True)
    parser.add_argument('-tf', '--target_embedding_file', type=str, required=True)
    parser.add_argument('-osf', '--output_source_file', type=str, required=True)
    parser.add_argument('-otf', '--output_target_file', type=str, required=True)
    parser.add_argument('-n', '--normalization', type=str, required=True,
                        choices=NORM_MODES)
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
    landmarks_B = obj_B.get('landmarks')
    weight_B = obj_B.get('weight')
    logger.info('- Target embedding: %s nodes' % len(nodes_A))

    assert( X_A.shape[1] == X_B.shape[1] )

    # Center the source/target, and color source with target covariance
    if args.normalization == CORAL_SRC_NORM:
        logger.info('[Aligning feature spaces with CORAL: source into target]')
        X_A, X_B = coral(X_A, X_B)

    # Center source/target, and color target non-landmarks with target landmarks covariance
    elif args.normalization == CORAL_TGT_NORM:
        logger.info('[Aligning feature spaces with CORAL: target non-landmarks into target landmarks]')

        # Extract the landmark/non-landmark indices
        node_to_index_B = dict(zip(nodes_B, range(len(nodes_B))))
        B_landmark_indices = [ node_to_index_B[n] for n in landmarks_B ]
        B_non_landmark_indices = [ node_to_index_B[n] for n in set(nodes_B)-set(landmarks_B) ]
        
        X_B_NL = X_B[B_non_landmark_indices]
        X_B_L = X_B[B_landmark_indices]
        
        # Use CORAL to align the non-landmarks to the landmarks, and
        # then reconstruct the full target space
        from sklearn.preprocessing import scale
        X_B_NL, X_B_L = coral(X_B_NL, X_B_L)
        X_B[B_landmark_indices] = X_B_L
        X_B[B_non_landmark_indices] = X_B_NL
        X_A = scale(X_A, with_std=False) # center source so it is comparable to the aligned target

    # No normalization
    elif args.normalization == NONE_NORM:
        logger.info('[Peforming no normalization]')
    else:
        raise NotImplementedError('Normalization mode "%s" not implemented.' % args.normalization)

    # Output aligned source embedding to file
    from util import save_embeddings
    save_embeddings(X_A, nodes_A, landmarks_A, args.output_source_file, weight_A)
    save_embeddings(X_B, nodes_B, landmarks_B, args.output_target_file, weight_B)
    
