

# Synetic Lethal (SL) interaction prediction with HANDL embeddings

Implementation of experiments in section 2.5 in [1]. In this experiment, we use HANDL embeddings to predict SL interactions in two different species of yeast and show that SLs are **co-located across species**. (Please see [1] for more details.)

From [1]:

> We test whether SLs are separated from non-SLs in Handl-space, and whether this separation extends across species, by training a classifier for gene pairs using Handl-embeddings as features. More specifically, we train a random forest (RF) to classify gene pairs as SLs or non-SLs within both species simultaneously, using the source embedding and the target embedded into source space. Without the Handl-embeddings, we would not be able to train a classifier for multiple species, since genes in different species would be in different vector spaces and have different dimensions. We perform 4-fold cross-validation, fixing the relative fraction of pairs from each species, and assess the degree of separation between SLs and non-SLs in Handl-space by evaluating the RF classifications with maximum F1 score (the harmonic mean of precision and recall), the area under the ROC curve (AUROC), and the area under the precision-recall curve (AUPRC). We report the average across the four folds, separating the results by species. We use a nested cross-validation strategy to choose the number of trees for the RF that maximizes the held-out AUPRC. For simplicity, all of our experiments in this section use S.c. as the source and S.p. as the target.


## Usage:

To obtain results for an experiment configured by a config file (as described below), simply execute:

	snakemake all --configfile <config file>

with a provided configuration file.

## Provided config files

We provide configuration files in `configs/` that correspond to experiments and results reported in [1].

*  `collins-roguev-rf.yml`, `collins-roguev-svm.yml`  - SL prediction using GI/SL data from Collins et. Al [2] and Roguev et. Al [3] for _S.c_ and _S.p_, respectively; using random forests and SVMs, respectively.
*  `biogrid.v3.4.157-rf.yml`, `biogrid.v3.4.157-svm.yml`-SL prediction using GI/SL data from BioGRID release v3.4.157 for _S.c_ and _S.p_, respectively; using random forests and SVMs, respectively. Non SLs are sampled from cartesian product of nodes obtained from PPI networks.
* `toy-rf.yml`, `toy-svm.yml`-SL prediction using synthetic data using random forests and SVMs, respectively. 
Note that inconclusive GIs are excluded
* `tiny.yml`- SL prediction using synthetic data using random forests with a single tree. Because the experiments in this directory can be time consuming, we include a configuration file to configure a run that executes quickly; this configuration can be used for sanity checking (and is used for continuous integration).

## Configuration

This experiment is implemented with [Snakemake](http://snakemake.readthedocs.io/en/stable/) and can be configured by a given [config file
](http://snakemake.readthedocs.io/en/stable/). Below is an example configuration file in JSON format (note that the comments included below are for annotation's sake and that a correctly written config file in JSON format cannot contain comments)

```javascript
{
"A": {
	"name": "<Name of species A>",
	"ppi": "<Path to species A PPI network>",
	"gi": "<Path to species A genetic interactions/SLs",
	"homologs": {
		"B": "<Path to homolog list (A to B)>"
		}
	},
"B": {
"name": "<Name of species B>",
"ppi": "<Path to species B PPI network>",
"gi": "<Path to species B genetic interations/SLs",
"homologs": {
	"A": "<Path to homolog list (B to A)"
	}
},
"n_landmarks": 400, // Number of landmarks for HANDL
"dataset_name": "<name of dataset (note that this is used as a prefix for the output directory)>",

// List of different cross-validation (C.V) strategies to test/run.
//     "pairs" --> C.V on gene pairs
//     "genes" --> C.V over individual genes
"hold_outs": ["genes", "pairs"],

// Different strategies for training classifier to test/run.
//     "both" --> train on both species simultaneously
//     "src"  --> train on source species (defaults to species A) only
//     "src"  --> train on target species (defaults to species B) only
"train_data": ["both", "src", "tgt"],

// Whether or not to sample negatives at random from pairs not seen in given GI files.
"sinatra_featurize": true,

// Number of cores/threads to use
"n_jobs": 42,

// Classifier to test/run. One of ["rf", "svm"]
"classifier": "rf"
}
```

### Provided config files

We provide configuration files in `configs/` that correspond to experiments and results reported in [1].

*  `collins-roguev.yml` - Configuration for SL prediction baseline using GI/SL data from Collins et. Al [2] and Roguev et. Al  [3] for _S.c_ and _S.p_ respectively.
*  `biogrid.v3.4.157`- Configuration for SL prediction baseline using GI/SL data from BioGRID [4] release v3.4.157 for _S.c_ and _S.p_ respectively. Non SLs are sampled from cartesian product of nodes obtained from PPI networks

Note that inconclusive GIs are excluded (Please see [1] section 4.4 for more details.)

## Reference

[1] Jason Fan, Anthony Cannistra, Inbar Fried, Tim Lim, Thomas Schaffner, Mark Crovella, Benjamin Hescott*, Mark DM Leiserson*. "A Multi-Species Functional Embedding Integrating Sequence and Network Structure."  _RECOMB 2018_(to appear)  [[bioRxiv]](https://www.biorxiv.org/content/early/2018/03/30/229211)  * equal contribution.

[2] Roguev, Sourav Bandyopadhyay, Martin Zofall, Ke Zhang, Tamas Fischer, et al. "Conservation and rewiring of functional modules revealed by an epistasis map in fission yeast".  _Science_, 322(5900):405–410, 8 2008.

[3] Sean R Collins, Kyle M Miller, Nancy L Maas, Assen Roguev, Jeffrey Fillingham, et al. "Functional dissection of protein complexes involved in yeast chromosome biology using a genetic interaction map".  _Nature_, 446(7137):806– 810, 7 2007

[4] Andrew Chatraryamontri, Rose Oughtred, Lorrie Boucher, Jennifer Rust, Christie Chang, et al. "The biogrid interaction database: 2017 update". _Nucleic Acids Research_, 45(D1):D369–D379, 2017.s