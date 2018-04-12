
# Plots to show close relationship between HANDL homology scores and cross-species functional similarity

In this example, we show that HANDL scores are correlated with functional similarity by producing  Figures 2a and 2b and/or Figures in Section S5 in [1].

## Usage

To produce Figures 2a and 2b or Figures in Section S5 in the Supplemental information from [1], simply execute:

	snakemake all --configfile <config file>

with a provided config file in `configs/`.

## Resnik scores

We have pre-computed resnik scores between several species and these are included in the default data set. Please see the [README](https://github.com/lrgr/HANDL) for more details.

## Provided configuration files:

We provide example configuration files for the following:
* `configs/sc-sp.yml`: a config file to produce plots of HANDL dissimilarity scores between _S.c_ and _S.p_ proteins projected into HANDL space with _S.c_ as the source and _S.p_ as the target. (See Section S5 in the Supplemental information from [1].)
* `configs/mouse-sc.yml` and `configs/human-sc/yml`: config files to produce plots between mouse and _S.c_ and human and _S.c_ in Section S5 in the Supplemental information from [1].
* `configs/human-mouse.yml`: config file to produce Figures 2a and 2b from [1].

## Scripts:
Python scripts to produce plots:
* [`plot_dissims.py`](https://github.com/lrgr/HANDL/blob/master/experiments/resnik-and-dissim-plots/plot_dissims.py): Script to plot HANDL dissimilarity scores of homolog pairs versus other pairs (see Figure 2b in [1]).
* [`plot_resnik_v_handl.py`](https://github.com/lrgr/HANDL/blob/master/experiments/resnik-and-dissim-plots/plot_resnik_v_handl.py): Script to plot pairs ranked by HANDL homology scores of homolog pairs versus Resnik scores (see Figure 2a in [1]).

## References
[1] Jason Fan, Anthony Cannistra, Inbar Fried, Tim Lim, Thomas Schaffner, Mark Crovella, Benjamin Hescott*, Mark DM Leiserson*. "A Multi-Species Functional Embedding Integrating Sequence and Network Structure."  _RECOMB 2018_(to appear)  [[bioRxiv]](https://www.biorxiv.org/content/early/2018/03/30/229211)  * equal contribution.