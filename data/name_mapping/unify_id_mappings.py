#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json
from uniprot_mapper import Uniprot_Mapper

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_files', type=str, required=True,
                    nargs='*', help='Uniprot ID mapping files (see https://goo.gl/uVgcwq)')
parser.add_argument('-o', '--output_file', type=str, required=True)
args = parser.parse_args( sys.argv[1:] )

# Load the files and output a unified version
mapper = Uniprot_Mapper(args.input_files)

# Output a unified map to file
output = dict(refseq=mapper.refseq_map,
              gene_ordered_locus=mapper.gene_ordered_locus_map,
              ensembl=mapper.ensembl_map,
              params=vars(args))
with open(args.output_file, 'w') as OUT:
    json.dump( output, OUT )
