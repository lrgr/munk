#!/usr/bin/env python

# Load required modules
import numpy as np
import gzip
import itertools
import argparse, os, sys, time

# Load our modules
from uniprot_mapper import Uniprot_Mapper

# the format of the file we are using for name mapping is at
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
# and the file is obtained from
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

# Parse command-line arguments
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--input_file', type=str, required=True)
    parser.add_argument('-of', '--output_file', type=str, required=True)
    parser.add_argument('-uf', '--uniprot_db_file', type=str, required=True)
    return parser
    
# for use in development 
def setup_args(in_f='tmpfile_homologs',
               out_f='foo.txt',
               uf = 'idmapping.dat.gz'):
    args = argparse.Namespace()
    args.input_file = in_f
    args.output_file = out_f 
    args.uniprot_db_file = uf
    return args

def run(args):
    # Set up logger
    sys.path.append(os.path.normpath(os.path.dirname(__file__) + '/../../src'))
    from logging_utils import getLogger
    logger = getLogger()

    # Load mapping and map
    logger.info('Loading uniprot mapping db')
    mapper = Uniprot_Mapper(args.uniprot_db_file, cached=True)
    n_in = 0
    n_out = 0

    logger.info('Starting mapping')
    with open(args.input_file, 'r') as in_f, open(args.output_file, 'w') as out_f:
        line = in_f.readline()
        while line:
            sc_gene, sp_gene = line.rstrip().split()
            n_in += 1

            sc_uni_name = mapper.map_refseq(sc_gene)
            sp_uni_name = mapper.map_refseq(sp_gene)

            if (sc_uni_name is not None) and (sp_uni_name is not None):
                n_out += 1
                out_f.write('{} {}\n'.format(sc_uni_name, sp_uni_name))
            line = in_f.readline()
    logger.info('Done. Read {} records.  Successfully mapped {} records.'.format(n_in, n_out))

if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
