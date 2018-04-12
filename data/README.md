# Default Data files for examples and experiments

Data used for examples and experiments in this project/repository have been pre-processed and hosted by our research group. The data required for all the experiments in this project can be downloaded with:

	snakemake all --configfile data.yml

Note that we include, as part of the default dataset,  pre-computed resnik scores that are large in size. If you do not need the resnik scores and would like to save disk space, you can instead run:
	
	snakemake all --configfile data-no-resnik.yml

Please see Section 4.4 of,

> Jason Fan, Anthony Cannistra, Inbar Fried, Tim Lim, Thomas Schaffner, Mark Crovella, Benjamin Hescott*, Mark DM Leiserson*. "A Multi-Species Functional Embedding Integrating Sequence and Network Structure."  _RECOMB 2018_(to appear)  [[bioRxiv]](https://www.biorxiv.org/content/early/2018/03/30/229211)  * equal contribution.

for complete details.


## A note on UniProt assession IDs

** Please note that for _S.c_ and _S.p_ some downloaded files have been pre-processed so that gene/protein names are mapped to UniProt assecssion IDs [1]. These files contain `-std` in their file names.**

## PPI networks for _S.p_, _S.c_, human and mouse.

PPI networks for _S.p_, _S.c_, human and mouse are downloaded to a directory named `ppi/`.
* `ppi/biogrid`: will contain ppi networks for _S.p_ and _S.c_ from the BioGRID database [1].
* `ppi/string-v9.1`: will contain ppi networks for human and mouse from the STRING database (v9.1) [2].

## Genetic interaction (GI) data for _S.p_, _S.c_

GI data will be downloaded to the following directories.
* `gi/biogrid`:  for BioGRID GI data for _S.p_ and _S.c_ [1].
* `gi/collins`: for GI data for _S.c_ from Collins et. Al [3].
* `gi/roguev`:  for GI data for _S.c_ from Roguev et. Al [4].

##  Homolog data/lists

Homolog lists (as processed from Homologene [5]) between the species, _S.c_, _S.p_, human, and mouse, will be downloaded to `homologs/`.

## Synthetic data

Synthetic (toy) data for two synthetic species will be downloaded to `toy/degree_sum_model/examples`. The synthetic PPI data is generated with the barabasi-albert model and two genes are determined to be SL if the sum of their degrees is greater than 14. The respective synthetic PPI networks and GI data can be found in:
* `toy/degree_sum_model/examples/alice_n100_d14`: for a species named _alice_ with 100 genes/proteins
* `toy/degree_sum_model/examples/bob_n120_d14`: for a species named _bob_ with 120 genes/proteins

## Pre-computed Resnik scores

Pre-computed Resnik scores between human-mouse, human-sc, and mouse-sc genes are downloaded to `resnik_scores`. Note that these files are large in size.

## References
[1] Andrew Chatraryamontri, Rose Oughtred, Lorrie Boucher, Jennifer Rust, Christie Chang, et al. "The biogrid interaction database: 2017 update". _Nucleic Acids Research_, 45(D1):D369–D379,

[2] Andrea Franceschini, Damian Szklarczyk, Sune Frankild, Michael Kuhn, Milan Simonovic, Alexander Roth, Jianyi Lin, Pablo Minguez, Peer Bork, Christian von Mering, and Lars J. Jensen. "String v9.1: protein-protein interaction networks, with increased coverage and integration". _Nucleic Acids Research_, 41(D1):D808–D815, 2013.

[3] Sean R Collins, Kyle M Miller, Nancy L Maas, Assen Roguev, Jeffrey Fillingham, et al. "Functional dissection of protein complexes involved in yeast chromosome biology using a genetic interaction map".  _Nature_, 446(7137):806– 810, 7 2007

[4] Assen Roguev, Sourav Bandyopadhyay, Martin Zofall, Ke Zhang, Tamas Fischer, et al. "Conservation and rewiring of functional modules revealed by an epistasis map in fission yeast".  _Science_, 322(5900):405–410, 8 2008.

[5] Eric W Sayers, Tanya Barrett, Dennis A Benson, Evan Bolton, Stephen H Bryant, et al. "Database resources of the national center for biotechnology information". _Nucleic acids research_, 39(suppl 1):D38–D51, 2011.