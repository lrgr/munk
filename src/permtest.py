import argparse

import networkx as nx
import numpy as np

from sklearn.externals.joblib import Parallel, delayed
from handl import separate_scores

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--source_ppi_edgelist',
                        type=str, required=True)
    parser.add_argument('-t', '--target_ppi_edgelist',
                        type=str, required=True)
    parser.add_argument('-sf', '--source_data_files', nargs='*', 
                        type=str, required=True)
    parser.add_argument('-tf', '--target_data_files', nargs='*', 
                        type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=str)
    parser.add_argument('-j', '--n_jobs', type=int, 
                        required=False, default=2)
    return parser.parse_args()

def difference_of_means(source_data_fp, target_data_fp):
    return 1.0

def main(args):
    source_data_files = args.source_data_files
    target_data_files = args.target_data_files
    assert(len(source_data_files) == len(target_data_files))

    print(source_data_files)
    print(target_data_files)
    print(args.source_ppi_edgelist)
    print(args.target_ppi_edgelist)
    print(args.n_jobs)

    n_permutations = len(source_data_files)
    source_target_data_files = zip(source_data_files, target_data_files)

    diffs = Parallel(n_jobs=args.n_jobs)(
                delayed(difference_of_means)(source_data_fp, target_data_fp)
                for source_data_fp, target_data_fp in 
                    source_target_data_files)
    print(diffs)

if __name__ == '__main__':
    main(parse_args())
