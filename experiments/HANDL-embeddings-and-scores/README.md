# HANDL Embeddings and HANDL homology scores

In this directory you will find example usage of [`/scripts/compute_embeddings.py`](https://github.com/lrgr/HANDL/blob/master/scripts/compute_embeddings.py) provided by the release of [HANDL](https://github.com/lrgr/HANDL). In this example we use Snakemake to configure the execution of `compute_embeddings.py`.

## Usage:

To compute HANDL embeddings and HANDL homology scores between _S.p_ and _S.c_ (fission and baker's yeast), simply execute:
	
	snakemake all --configfile configs/sc-sp.yml

The following files will be outputted in `output/sc-to-sp/`for HANDL embeddings form _S.c_ (source species) to _S.p_ (target species):
* `sc-to-sp-HANDL-landmarks.tsv` : a tab separated list of landmark pairs (homolog pairs between species) used to compute HANDL embeddings
* `sc-to-sp-HANDL-similarity-scores.pkl`: a pickled dictionary containing HANDL homology scores between _S.c_ and _S.p_
* `sc-to-sp-runtimes.json`: runtimes for computing HANDL embeddings and HANDL homology scores
* `sc-to-sp_sc-HANDL-embedding.pkl`: a pickled dictionary containing HANDL embeddings for _S.c_
* `sc-to-sp_sp-HANDL-embedding.pkl`: a pickled dictionary containing HANDL embeddings for _S.p_

Files following similar naming schemes containing HANDL Embeddings and HANDL homology scores with _S.p_ as the source and _S.c_ as the target species are also outputted to `output/sp-to-sc`. Note that outputted, pickled dictionaries are saved, and can be loaded, using Scikit-learn's [joblib](http://pythonhosted.org/joblib/) module.

### Configuration:

Execution of the experiment/script is implemented/written with [Snakemake](http://snakemake.readthedocs.io/en/latest/index.html) and is configured by a given [config file] http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html?highlight=config%20file.

The provided configuration file should provide information for a pair of species (A and B), and the included Snakefile will compute HANDL embeddings from A to B, and from B to A.

A configuration file can be defined as follows:
```YAML
# An example configuration file in YAML format for HANDL embeddings between S.c and S.p

# Name of species A
species-A: sc

# Name of species B
species-B: sp

# Path to PPI edgelist corresponding to species A
A-ppi: ../../data/ppi/biogrid/sc/sc-biogrid.v3.4.157-ppi-std.tsv

# Path to PPI edgelist corresponding to species B
B-ppi: ../../data/ppi/biogrid/sp/sp-biogrid.v3.4.157-ppi-std.tsv

# Path to tab separated homolog lists for HANDL embeddings from A (source) to B (target)
A-to-B-homologs: ../../data/homologs/sc-sp/sc-sp-homologs.txt

# Path to tab separated homolog lists for HANDL embeddings from B (source) to A (target)
B-to-A-homologs: ../../data/homologs/sc-sp/sp-sc-homologs.txt

# Number of landmarks to use for HANDL embeddings
n_landmarks: 400
```

