
# Synetic Lethal (SL) interaction prediction baseline with homologs

In this experiment, we produce a baseline for SL prediction using homologs. In particular, we predict homologs with the following scheme. We say that for a pair of genes (a,b) in a source species, and a pair of genes (u,v) in a target species, that (u,v) is SL iff (a,b) is SL and (a,u) and (b,v) are homologs.

Please see [1] section 2.5 for more details.



## Usage:
To obtain results for an experiment configured by a config file (as described below), simply execute:

	snakemake --configfile <config file>

## Configuration

This experiment is implemented with [Snakemake](http://snakemake.readthedocs.io/en/stable/) and can be configured by a given [config file
](http://snakemake.readthedocs.io/en/stable/). Below is an example configuration file in YAML format.

```YAML
# Example configuration file to configure SL prediction baseline experiment with data 
# from Collins et. Al. and Roguev et Al. for SL predictions between S.c and S.p, respectively.
# SL predictions are transferred between and evaluated for given species A and B
#

# Run-name - Note that this is used as a prefix for the results output file.
run-name: biogrid.v3.4.157

# Species A information
A:
  # Data source for A (used in outputted results table)
  data-name: biogrid.v3.4.157
  # Species A name
  name: sc
  # Path to species A genetic interactions (GI) file 
  # (See /data/ for more information for GI file formats)
  gi: ../../data/gi/biogrid/sc/sc-biogrid.v3.4.157-sls-std.tsv
  # Path to species B PPI network edglist
  ppi: ../../data/ppi/biogrid/sc/sc-biogrid.v3.4.157-ppi-std.tsv

B:
  data-name: biogrid.v3.4.157
  name: sp
  gi: ../../data/gi/biogrid/sp/sp-biogrid.v3.4.157-sls-std.tsv
  ppi: ../../data/ppi/biogrid/sp/sp-biogrid.v3.4.157-ppi-std.tsv

# Path to A to B tab separated list of homolog pairs
homologs: ../../data/homologs/sc-sp/sc-sp-homologs.txt

# Flag for whether to sample non-SL interactions at random from pairs found in PPI networks
# that do not have measured SL interactions.
# If True, non-SLs are sampled. If false, *measured* non-SLs are obtained from given GI files
sample-negs: True
#sample-negs: False
```

### Provided config files

We provide configuration files in `configs/` that correspond to experiments and results reported in [1].

*  `collins-roguev.yml` - Configuration for SL prediction baseline using GI/SL data from Collins et. Al [2] and Roguev et. Al [3] for _S.c_ and _S.p_ respectively.
*  `biogrid.v3.4.157`- Configuration for SL prediction baseline using GI/SL data from BioGRID [4] release v3.4.157 for _S.c_ and _S.p_ respectively. Non SLs are sampled from cartesian product of nodes obtained from PPI networks

Note that inconclusive GIs are excluded (Please see [1] section 4.4 for more details.)

## Reference

[1] Jason Fan, Anthony Cannistra, Inbar Fried, Tim Lim, Thomas Schaffner, Mark Crovella, Benjamin Hescott*, Mark DM Leiserson*. "A Multi-Species Functional Embedding Integrating Sequence and Network Structure."  _RECOMB 2018_(to appear)  [[bioRxiv]](https://www.biorxiv.org/content/early/2018/03/30/229211)  * equal contribution.

[2] Roguev, Sourav Bandyopadhyay, Martin Zofall, Ke Zhang, Tamas Fischer, et al. "Conservation and rewiring of functional modules revealed by an epistasis map in fission yeast".  _Science_, 322(5900):405–410, 8 2008.

[3] Sean R Collins, Kyle M Miller, Nancy L Maas, Assen Roguev, Jeffrey Fillingham, et al. "Functional dissection of protein complexes involved in yeast chromosome biology using a genetic interaction map".  _Nature_, 446(7137):806– 810, 7 2007

[4] Andrew Chatraryamontri, Rose Oughtred, Lorrie Boucher, Jennifer Rust, Christie Chang, et al. "The biogrid interaction database: 2017 update". _Nucleic Acids Research_, 45(D1):D369–D379, 2017.s