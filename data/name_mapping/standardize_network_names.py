#!/usr/bin/env python

# Load required modules
import numpy as np
import gzip
import itertools
import argparse, os, sys, time

# Load our modules and set up logger
from uniprot_mapper import Uniprot_Mapper

# the format of the file we are using for name mapping is at
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
# and the file is obtained from
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

# Parse command-line arguments
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--input_file', type=str, required=True, action='append')
    parser.add_argument('-of', '--output_file', type=str, required=True, action='append')
    parser.add_argument('-uf', '--uniprot_db_file', type=str, required=True)
    parser.add_argument('-sp', '--species', type=str, required=True)
    return parser

def announce(message):
    print(time.strftime('%H:%M:%S'),message)
    sys.stdout.flush()

def map_file(input_file, output_file, err_file, map_name):

    n_in = 0
    n_out = 0
    with open(input_file, 'r') as in_f, open(output_file, 'w') as out_f, open(err_file, 'w') as err_f:
        line = in_f.readline()
        while line:
            fields = line.rstrip().split()
            left_gene, right_gene = fields[:2]
            n_in += 1

            left_uni_name = map_name(left_gene)
            right_uni_name = map_name(right_gene)

            if (left_uni_name is not None) and (right_uni_name is not None):
                n_out += 1
                out_f.write('{} {} {}\n'.format(left_uni_name,
                                                right_uni_name,
                                                ' '.join(fields[2:])))
            else:
                if left_uni_name is None:
                    err_f.write('{} '.format(left_gene))
                if right_uni_name is None:
                    err_f.write('{} '.format(right_gene))
                err_f.write('\n')
            line = in_f.readline()
    return n_in, n_out

def run(args):
    # Set up logger
    sys.path.append(os.path.normpath(os.path.dirname(__file__) + '/../../src'))
    from logging_utils import getLogger
    logger = getLogger()

    # Load mapping and map
    logger.info('Loading uniprot mapping db')
    mapper = Uniprot_Mapper(args.uniprot_db_file, cached=True)

    if args.species == 'sc':
        map_name = mapper.map_gene_ordered_locus
    elif args.species == 'sp':
        map_name = mapper.map_ensembl

    for in_f, out_f in zip(args.input_file, args.output_file):
        logger.info('Starting mapping {} -> {}'.format(in_f, out_f))
        err_f = in_f + '.errs'
        n_in, n_out = map_file(in_f, out_f, err_f, map_name)
        logger.info('Done. Read {} records.  Successfully mapped {} records.'.format(
            n_in, n_out))

if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
