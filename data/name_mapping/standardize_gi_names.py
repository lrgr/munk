#!/usr/bin/env python

# Load required modules
import numpy as np
import gzip
import itertools
import argparse, os, sys, time

# Load our modules
from uniprot_mapper import Uniprot_Mapper
sys.path.append(os.path.normpath(os.path.dirname(__file__) + '/../gi'))
from gi_io import output_gi_dataframe

# the format of the file we are using for name mapping is at
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
# and the file is obtained from
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

# Parse command-line arguments
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--input_file', type=str, required=True, action='append')
    parser.add_argument('-of', '--output_file_prefix', type=str, required=True, action='append')
    parser.add_argument('-uf', '--uniprot_db_file', type=str, required=True)
    parser.add_argument('-sp', '--species', type=str, required=True)
    return parser

def map_file(input_file, output_file_prefix, map_name, test=False):
    import pandas as pd
    err_file = output_file_prefix + '.errs'
    n_in = 0
    n_out = 0

    # Quick and dirty entrypoint for testing
    if test:
        errors = [{'Gene': 'dummy'}]
    else:
        errors = []

    items = []
    df = pd.DataFrame.from_csv(input_file, sep='\t', header=0, index_col=None)
    for i, r in df.iterrows():
        geneA, geneB = r['Gene A'], r['Gene B']
        score, category= r['Score'], r['Category']
        uni_name_A = map_name(geneA)
        uni_name_B = map_name(geneB)
        if (uni_name_A is not None) and (uni_name_B is not None):
            items.append({
                 'Gene A': uni_name_A,
                 'Gene B': uni_name_B,
                 'Score': float(score),
                 'Category' : category})
        else:
            if uni_name_A is None:
                errors.append({'Gene': geneA})
            if uni_name_B is None:
                errors.append({'Gene': geneB})

    std_df = pd.DataFrame(items)[['Gene A', 'Gene B', 'Score', 'Category']]
    std_df.sort_values(['Gene A', 'Gene B'])

    # Save a grep-able TSV file and a Numpy file with gene index
    output_gi_dataframe(std_df, output_file_prefix, 'tsv')
    output_gi_dataframe(std_df, output_file_prefix, 'npy')

    if len(errors) > 0:

        err_df = pd.DataFrame(errors)[['Gene']]
        err_df.sort_values(['Gene'])
        err_df.to_csv(err_file, index=False, sep='\t')
        n_out = len(errors)

    n_in = df.shape[0]
    n_out = len(items)
    return n_in, n_out

def run(args):
    # Set up logger
    sys.path.append(os.path.normpath(os.path.dirname(__file__) + '/../../src'))
    from logging_utils import getLogger
    logger = getLogger()

    # Load mapping and map
    test = False
    map_name = None

    if args.species == 'test':
        test = True
        map_name = lambda n : 'std_name'
    else:
        mapper = Uniprot_Mapper(args.uniprot_db_file, cached=True)
        if args.species == 'sc':
            map_name = mapper.map_gene_ordered_locus
        elif args.species == 'sp':
            map_name = mapper.map_ensembl
        else:
            raise NotImplementedError()


    for in_f, out_f_prefix in zip(args.input_file, args.output_file_prefix):
        logger.info('Starting mapping {} -> {}'.format(in_f, out_f_prefix + '.tsv'))
        n_in, n_out = map_file(in_f, out_f_prefix, map_name, test=test)
        logger.info('Done. Read {} records.  Successfully mapped {} records.'.format(
            n_in, n_out))

if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
