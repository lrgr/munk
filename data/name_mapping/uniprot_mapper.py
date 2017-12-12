import numpy as np
import gzip
import itertools
import os, sys, json

# the format of the file we are using for name mapping is at
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
# and the file is obtained from
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

# SGD ORFS are Gene_OrderedLocusName names in uniprot
# SPO systematic names are EnsemblGenome names in uniprot
# Homologene 'HG-Gene_prot-acc' names are RefSeq names in uniprot

class Uniprot_Mapper:

    def __init__(self, file_names, cached=False):
        # Convert a single file name to a list of files
        if type(file_names) == type(''):
            file_names = [file_names]
                
        # Figure out if we're looking at a cache or not
        if cached:
            with open(file_names[0], 'r') as IN:
                data = json.load(IN)
                self.refseq_map = data['refseq']
                self.gene_ordered_locus_map = data['gene_ordered_locus']
                self.ensembl_map = data['ensembl']
        else:
            self.refseq_map = {}
            self.gene_ordered_locus_map = {}
            self.ensembl_map = {}
            
            for uniprot_db_file in file_names:
                with gzip.open(uniprot_db_file,'r') as in_f:
                    for line in in_f:
                        (uniprot_ref, identifier, value) = line.decode('utf-8').rstrip().split('\t')
                        if identifier == 'RefSeq':
                            self.refseq_map[value] = uniprot_ref
                        elif identifier == 'Gene_OrderedLocusName':
                            self.gene_ordered_locus_map[value] = uniprot_ref
                        elif identifier == 'EnsemblGenome':
                            self.ensembl_map[value] = uniprot_ref

    # dont want to use exceptions for signaling here
    # because they are slow and will be raised a lot
    # so caller should check if return value is not None
    @staticmethod
    def __lookup(name, mapper):
        if name in mapper:
            return mapper[name]
        # this special case addresses a coding issue in Pombe
        # labs tend to report Pombe systematic IDs with trailing capital 'C'
        # but correct encoding uses lowercase 'c'
        elif name[-1] == 'C':
            name = name[:-1] + 'c'
            if name in mapper:
                return mapper[name]
            else:
                return None
        else:
            return None
        
    def map_refseq(self, name):
        return self.__lookup(name, self.refseq_map)

    def map_gene_ordered_locus(self, name):
        return self.__lookup(name, self.gene_ordered_locus_map)

    def map_ensembl(self, name):
        return self.__lookup(name, self.ensembl_map)

