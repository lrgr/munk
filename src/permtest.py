import argparse
import json
import networkx as nx
import numpy as np
from sklearn.externals import joblib
from sklearn.externals.joblib import Parallel, delayed

import handl
import util
from i_o import get_logger

def parse_args():
    parser = argparse.ArgumentParser()

    # arguments to perfrom HANDL
    parser.add_argument('-s', '--source_edgelist',
                        type=str, required=True)
    parser.add_argument('-t', '--target_edgelist',
                        type=str, required=True)
    parser.add_argument('-hom', '--homolog_list', 
                        type=str, required=True)
    parser.add_argument('-nl', '--n_landmarks',
                        type=int, required=False)
    parser.add_argument('-lam', '--lam', 
                        type=float, required=False)

    parser.add_argument('-sf', '--source_data_files', nargs='*', 
                        type=str, required=True)
    parser.add_argument('-tf', '--target_data_files', nargs='*', 
                        type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-dof', '--diffs_output_file', type=str, required=True)
    parser.add_argument('-j', '--n_jobs', type=int, 
                        required=False, default=2)
    return parser.parse_args()

def difference_of_means(source_data_fp, target_data_fp):
    return 1.0

def handl_embed_wrapper(source_data, 
                        target_data, 
                        homologs, 
                        n_landmarks,
                        ref_source_nodes,
                        ref_target_nodes):  
    source_C = source_data['C']
    target_D = target_data['D']
    source_nodes = source_data['nodes']
    target_nodes = target_data['nodes']
    assert(source_nodes == ref_source_nodes)
    assert(target_nodes == ref_target_nodes)
    source_n2i = dict((n,i) for i, n in enumerate(source_nodes))
    target_n2i = dict((n,i) for i, n in enumerate(target_nodes))
    homolog_idxs =  [(source_n2i[s_n], target_n2i[t_n]) for
                     s_n, t_n in homologs]
    landmark_idxs = homolog_idxs[:n_landmarks]
    source_C, target_C_hat = handl.handl_embed_matrices(source_C, target_D, landmark_idxs)
    sim_scores = source_C @ target_C_hat.T
    return sim_scores, landmark_idxs, homolog_idxs

def dissim_diff(sim_scores, landmark_idxs, homolog_idxs):
    dissim_scores = 1. / sim_scores
    _,_,_, hom_hom_scores, other_scores = \
        handl.separate_scores(dissim_scores, landmark_idxs, homolog_idxs)
    other_mean = np.mean(other_scores)
    hom_mean = np.mean(hom_hom_scores)
    return other_mean - hom_mean, other_mean, hom_mean

def difference_in_means(source_data,
                        target_data,
                        homologs,
                        n_landmarks,
                        ref_source_nodes,
                        ref_target_nodes):
    handl_results = handl_embed_wrapper(source_data,
                                        target_data,
                                        homologs,
                                        n_landmarks,
                                        ref_source_nodes,
                                        ref_target_nodes)
    return dissim_diff(*handl_results)

def difference_in_means_from_files(source_data_file,
                        target_data_file,
                        homologs,
                        n_landmarks,
                        ref_source_nodes,
                        ref_target_nodes):
    try:
        source_data = joblib.load(source_data_file)
        target_data = joblib.load(target_data_file)
        return difference_in_means(source_data,
                                   target_data,
                                   homologs,
                                   n_landmarks,
                                   ref_source_nodes,
                                   ref_target_nodes)
    except:
        print('ERROR with {}, {}'.format(source_data_file, target_data_file))
        return None

def one_tail_pval(observed, dist):
    n_less_than = np.sum(np.where(dist > observed))
    return n_less_than / len(dist)

def effect_size(observed, control):
    std_dev = np.std(control)
    control_mean = np.mean(control)
    return (observed - control_mean) / std_dev
    
def main(args):
    log = get_logger()
    source_data_files = args.source_data_files
    target_data_files = args.target_data_files
    assert(len(source_data_files) == len(target_data_files))
    
    log.info('Loading homologs list from %s', args.homolog_list)
    raw_homologs = util.read_homolog_list(args.homolog_list)

    log.info('Loading source edgelist from %s', args.source_edgelist)
    source_G = nx.read_edgelist(args.source_edgelist, encoding='ascii')
    log.info('Loading target edgelist from %s', args.target_edgelist)
    target_G = nx.read_edgelist(args.target_edgelist, encoding='ascii')

    log.info('Computing HANDL embeddings real PPI networks with %d landmarks', args.n_landmarks)
    n_landmarks = args.n_landmarks

    source_G = util.simple_two_core(source_G)
    target_G = util.simple_two_core(target_G)
    homologs = handl.homologs_in_graphs(source_G, target_G, raw_homologs)

    source_nodes = sorted(source_G.nodes())
    target_nodes = sorted(target_G.nodes())

    source_D = handl.regularized_laplacian(source_G, source_nodes, args.lam)
    target_D = handl.regularized_laplacian(target_G, target_nodes, args.lam)
    source_C = handl.rkhs_factor(source_D)

    source_data = dict(C=source_C, nodes=source_nodes)
    target_data = dict(D=target_D, nodes=target_nodes)

    real_diff, other_mean, hom_mean = difference_in_means(source_data, 
                                    target_data, 
                                    homologs, 
                                    n_landmarks,
                                    source_nodes,
                                    target_nodes)
    log.info('Difference between homolog and non-homolog similarity scores %f', real_diff)

    n_permutations = len(source_data_files)
    log.info('Loading %d pairs of graphs from disk... and computing differences in means...', n_permutations)
    means = Parallel(n_jobs=args.n_jobs)(
                delayed(difference_in_means_from_files)
                       (source_data_fp, 
                        target_data_fp,
                        homologs,
                        n_landmarks,
                        source_nodes,
                        target_nodes)
                for source_data_fp, target_data_fp in 
                    zip(source_data_files, target_data_files))
    errs = [1 for m in means if m is None]
    print('ERRORS:', len(errs))
    means = [m for m in means if m is not None]

    rand_G_mean_diffs, rand_G_other_means, rand_G_hom_means = zip(*means)
    log.info('Mean difference in means between homologs and non-homologs scores for random graphs %f', np.mean(rand_G_mean_diffs))
    p_val = one_tail_pval(real_diff, rand_G_mean_diffs)
    e_size = effect_size(real_diff, rand_G_mean_diffs)

    log.info('P-value: %f', p_val)
    log.info('Effect-size: %f', e_size)
    log.info('# permutations: %d', n_permutations)
    results = dict(pval=p_val, effect_size=e_size, n_permutations=n_permutations, errs=len(errs))

    log.info('Writing results to %s', args.output_file)
    with open(args.output_file, 'w') as OUT:
        json.dump(results, OUT, indent=2)
    diffs = dict( real=dict(mean_diff=real_diff, 
                            non_hom_mean=other_mean, 
                            hom_mean=hom_mean),
                  random=dict(mean_diffs=rand_G_mean_diffs, 
                              non_hom_means=rand_G_other_means, 
                              hom_means=rand_G_hom_means))
    log.info('Writing values of differences to %s', args.diffs_output_file)
    joblib.dump(diffs, args.diffs_output_file)

if __name__ == '__main__':
    main(parse_args())
