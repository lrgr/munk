#!/usr/bin/env python

# Load required modules
import sys, os, argparse, networkx as nx

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-if', '--input_file', type=str, required=True)
parser.add_argument('-of', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

# Load the input file
G = nx.Graph()
with open(args.input_file, 'r') as IN:
    header = IN.readline().rstrip('\n').split('\t')
    for l in IN:
        arr = l.rstrip('\n').split('\t')
        p1 = arr[2].split()[-1].split(':')[1]
        p2 = arr[3].split()[-1].split(':')[1]
        G.add_edge(p1, p2)

# Find the largest connected component
ccs = sorted(list(nx.connected_components(G)), key=lambda cc: len(cc), reverse=True)
H = G.subgraph(ccs[0])
print('* Loaded graph with %s connected components...' % len(ccs))
print('\t- %s nodes (%s in largest CC)' % (G.number_of_nodes(), H.number_of_nodes()))
print('\t- %s edges (%s in largest CC)' % (G.number_of_edges(), H.number_of_edges()))

# Output to file
nx.write_edgelist(G, args.output_file, data=True, delimiter='\t')
