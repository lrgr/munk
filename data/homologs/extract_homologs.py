import numpy as np
from xml.etree import ElementTree as ET
import gzip
import itertools
import argparse, os, sys

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-if', '--input_file', type=str, default='homologene.xml.gz')
parser.add_argument('-of', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

# takes 3-4 mins and needs lots of memory
tree = ET.parse(gzip.open(args.input_file))
root_element = tree.getroot()

# Root contains one entry: HG-EntrySet_entries
assert(len(root_element) == 1)
entries = root_element[0]

# Parsing is guided by documentation in 
# ftp://ftp.ncbi.nih.gov/pub/HomoloGene/HomoloGene_Field_Description.txt

# from ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/build_inputs/taxid_taxname
Sp_taxid = '4896'
Sc_taxid = '4932'

n_pairs = 0

with open(args.output_file, 'w') as of:

    # each entry is one group of homologous genes
    # there are multiple species present and there
    # can be multiple genes for a give species in the group
    for entry in entries:

        # set of Gene objects 
        genes = entry.find('HG-Entry_genes')
        # set of Stats objects
        stats = entry.find('HG-Entry_distances')

        # what we want is the following.  if there is a group in which both
        # Sp and Sc genes are present, choose the Sp-Sc pair which is reciprocal-best
        # in BLAST score.  This is stored in the Stats_recip-best field
        # in the stats element

        Sp_genes = []
        Sc_genes = []

        for gene in genes:

            # Contains information about each gene in a HomoloGene group.  Contains
            # the following:
            # - Gene_geneid: GeneID (from EntrezGene) identifier
            # - Gene_symbol: "best" symbol for gene (official if available)
            # - Gene_title: name of the gene
            # - Gene_taxid: taxonomic id for the gene
            # - Gene_prot-gi: gi of protein used for calculations
            # - Gene_prot-acc: accession of protein used for calculations
            # - Gene_prot-len: length of protein used for calculations
            # - Gene_nuc-gi: gi of corresponding mRNA to protein used in
            #   calculations (if available)
            # - Gene_nuc-acc: accession of corresponding mRNA to proteinused in
            #   calculations (if available)
            # - Gene_domains: set of Domain objects

            taxid = gene.find('HG-Gene_taxid').text
            if taxid == Sp_taxid:
                Sp_genes.append(gene)
            elif taxid == Sc_taxid:
                Sc_genes.append(gene)

        if (Sp_genes and Sc_genes):

            # create a structure that allows us to quickly look up stats for a gene
            # pair regardless of ordering (XML representation is inflexible here)
            stats_array = {frozenset([stat.find('HG-Stats_gi1').text, 
                            stat.find('HG-Stats_gi2').text]): stat for stat in stats}

            for Sp, Sc in itertools.product(Sp_genes, Sc_genes):

                stat = stats_array[frozenset([Sc.find('HG-Gene_prot-gi').text,
                                              Sp.find('HG-Gene_prot-gi').text])]


                if stat.find('HG-Stats_recip-best').attrib['value'] == 'true':
                    of.write('{} {}\n'.format(                    
                       Sc.find('HG-Gene_prot-acc').text,
                       Sp.find('HG-Gene_prot-acc').text))
                    n_pairs += 1

                # stats object:
                # Contains different pairwise statistics between the proteins in a
                # HomoloGene group.  Contains the following:

                # - Stats_gi1: gi of the first protein in the pair
                # - Stats_gi2: gi of the second protein in the pair
                # - Stats_nuc-change: ratio of nucleotide differences between the pair
                # - Stats_nuc-change-jc: ratio of nucleotide differences between the
                #   pair (corrected for back substitions through Jukes and Cantor
                #   formula)
                # - Stats_prot-change: ratio of amino acid differences between the pair
                # - Stats_ka: Ka for the pair (ratio of non-synonymous differences per
                #   non-synonymous site)
                # - Stats_ks: Ks for the pair (ratio of synonymous differences per
                #   synonymous site)
                # - Stats_knr: Knr for the pair (ratio of radical non-synonymous
                #   differences per radical non-synonymous site)
                # - Stats_knc: Knc for the pair (ratio of conservative non-synonymous
                #   differences per conservative non-synonymous site)
                # - Stats_recip-best: is this pair reciprocal best using bitscore as
                #   the measurement
