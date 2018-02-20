#!/usr/bin/env python

################################################################################
# SETUP
################################################################################
# Load required modules
import sys, os, argparse, logging

# Load our modules
this_dir = os.path.dirname(__file__)
sys.path.append(os.path.join(this_dir, '../../src'))
from i_o import get_logger

################################################################################
# MAIN
################################################################################
# Parse command-line arguments
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ef', '--embedding_file', type=str, required=True)
    parser.add_argument('-gif', '--gi_file', type=str, required=True)
    parser.add_argument('-of', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbosity', type=int, required=False,
                        default=logging.INFO)
    return parser

# Main
def run( args ):
    # Create logger
    logger = get_logger(args.verbosity)

    # Output to file
    with open(args.output_file, 'w') as OUT:
        OUT.write('done!\n')

if __name__ == '__main__':
    run( get_parser().parse_args(sys.argv[1:]) )
