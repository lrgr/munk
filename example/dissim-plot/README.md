# Homologs are more similar in HANDL space than non-homologs

In this example, we show that HANDL scores are correlated with functional similarity by reproducing Figure 2 from [REF].

### Usage

To produce Figure 2 from [REF] simply execute `snakemake all --configfile configs/<configfile>`

### Provided configuration files:

We provide example configuration files for the following:
* `configs/sc-sp` is a configuration file to produce plots of HANDL dissimilarity scores between _S.c_ and _S.p_ proteins projected into HANDL space with _S.c_ as the source and _S.p_ as the target, respectively. PPI networks from which HANDL scores and embeddings are computed can be found in the `data` directory in the root folder and are obtained from BioGrid v3.4.157.
