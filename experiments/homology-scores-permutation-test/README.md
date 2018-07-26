# Permutation tests showing that MUNK captures shared biological information beyond node degree

See [1] section 2.3 for more complete details.

This experiment/directory contains implementation and scripts used to perform permutation tests in [1] section 2.3.

The experiment is described in [1] as follows:
> "We assess the statistical significance of the difference in MUNK
> homology scores between homologous and non-homologous pairs by
> generating 1000 pairs of random networks in which each node is given
> degree very close to that in the original network, but in which edges
> have been randomized. Specifically, we follow the method of Newman et.
> al [2] to generate graphs with given degree distributions and remove
> self loops and parallel edges afterwards. We then compute an empirical
> P-value by counting the number of pairs of random networks for which
> the difference in mean MUNK homology scores between homologs and
> non-homologs is greater than that observed in real PPI networks."

## Usage

This permutation test is configured and implemented with [Snakemake](http://snakemake.readthedocs.io/en/stable/). To run this experiment, simply execute:

	snakemake all --configfile configs/sc-sp.yml

We also include the option to speed up the process of generating networks for our permutation test by including a scipt to generate permuted networks (quickly) in parallel using a SLURM cluster. If you would like to generate permuted networks quickly and in parallel before you run the permutation test with Snakemake (as described above), you first can execute the following command **before** you run `snakemake all --configfile configs/sc-sp.yml`:

	sh gen_graphs_with_slurm.sh configs/sc-sp.yml

Morevoer, because generating 1000 random networks, as in [1], can be time consuming, we have included other config files that configure permutation test runs for sanity checking/testing. These runs can be performed with synthetic data (using `configs/alice-bob.yml`) or with fewer permuted networks (with `configs/sc-sp-test.yml`).

## Scripts:

Python scripts:
* [`gen_random_network.py`](https://github.com/lrgr/munk/blob/master/experiments/permutation-test/gen_random_network.py): A script to generate and save a single random network conditioned on the node degree distribution of a given "real" network. The parameters are as follows:
	* `-e`, `--edgelist`: PPI edgelist of seed/real network to permute/randomize
	* `-o`, `--output`: Output path
	* `-n`, `--n_tries`: Number of attempts to try to generate random network. (Note that sometimes the random network generation procedure will generate networks with fewer nodes than that of the the real/seed network. We count these instances as errors and this parameter sets how many times the script should repeat until it successfully generates a random/permuted network.)
* [`permtest.py`](https://github.com/lrgr/munk/blob/master/experiments/permutation-test/gen_random_network.py): A script to run the permutation test.
	* `-s'`, `--source_edgelist`: PPI edgelist of the "real" source species PPI network
	* `-t'`, `--target_edgelist`: PPI edgelist of the "real" target species PPI network
	* `-hom'`, `--homolog_list`: Tab separated list of homolog pairs between source and target species
	* `-nl'`, `--n_landmarks`: Number of landmarks for computing MUNK homology scores
	* `-lam`, `--lam`: $\lambda$ value to use for regularized laplacians within MUNK
	* `-sf'`, `--source_data_files`: list of randomized/permuted source networks (as outputted by [`gen_random_network.py`](https://github.com/lrgr/munk/blob/master/experiments/permutation-test/gen_random_network.py)).
	* `-tg'`, `--target_data_files`: list of randomized/permuted target networks (as outputted by [`gen_random_network.py`](https://github.com/lrgr/munk/blob/master/experiments/permutation-test/gen_random_network.py)).
	* `-o'`, `--output_file`: path to save results.
	* `-dof'`, `--diffs_output_file`: path to pickle detailed results/data
	* `-j'`, `--n_jobs`: number of cores/threads to use

Shell script for generating random networks using a SLURM cluster
* [`gen_graphs_with_slurm.sh`](https://github.com/lrgr/munk/blob/munk-package/experiments/permutation-test/gen_graphs_with_slurm.sh): Use a cluster, as defined in `cluster.yml` to generate many random networks in parallel!
	* To generate random networks for a permutation test configured with a configuration file, simply execute: `gen_graphs_with_slurm.sh <config file>`.

## References
[1] Jason Fan, Anthony Cannistra, Inbar Fried, Tim Lim, Thomas Schaffner, Mark Crovella, Benjamin Hescott*, Mark DM Leiserson*. "A Multi-Species Functional Embedding Integrating Sequence and Network Structure."  _RECOMB 2018_(to appear)  [[bioRxiv]](https://www.biorxiv.org/content/early/2018/03/30/229211)  * equal contribution.
[2] M. E. J. Newman, S. H. Strogatz, and D. J. Watts. "Random graphs with arbitrary degree distributions and their applications". _Phys. Rev. E_, 64:026118, Jul 2001.
