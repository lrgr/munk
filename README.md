# Homology Assessment across Networks using Diffusion and Landmarks (HANDL)

### HANDL
HANDL is an algorithm for embedding proteins in a target network (e.g. mouse) into a source network (e.g. from human). The HANDL algorithm was developed by Mark Crovella (Boston University), Benjamin Hescott (Northeastern University), and Max Leiserson (University of Maryland, College Park), and their respective research groups.

-----

### Setup
#### Dependencies
**Required Programs:** GNU Make, wget, Python 3.4.2+, PIP, conda

**Python Packages:** Install required modules with `conda env create -f environment.yml`

**Conda Environment:** Run `source activate HANDL` to work in `HANDL` conda environment

#### Data
HANDL requires as input a source and target PPI network, and a list of homologs mapping a subset of the nodes in each network. We include scripts for downloading and processing data for _S. cerevisiae_ (_Sc_) and _S. pombe_ (_Sp_) networks and homologs, and mapping them into the same namespace.

You can process the _Sc_ and _Sp_ data with the following commands.

1. **Download and process UniProt accession ID mappings.** Run `make all` in `data/name_mapping`.

2. **Download and process Sc and Sp homologs from Homlogene.** Run `make all` in `data/homologs`.

3. **Download and process BioGRID PPI networks.** Run `make all` in `data/ppi/biogrid`.

-----
### Usage - Scripts and command-line arguments
**Compute RKHS factors/factorization of regularized Laplacian of PPI networks with:** `factorized_laplacian.py`

This script computes and saves the regularized Laplacian and factored RKHS embedding of a given PPI network. 


Required parameters:

*   `-e`, `--edge_file` : Path to PPI network edge list
*   `-o`, `-output_file` : Path to where RKHS embedding should be saved
*   `-df`, `--diffusion_file` : Path to where regularized Laplacian should be saved
*   `-l`, `--lam` : Value of lambda in with respect to the regularized Laplacian

Outputs:
The regularized Laplacian is saved to a Python dictionary with the following keys and values:
* `D` :  A NumPy array containing the values of the regularized Laplacian
* `nodes`: List of nodes/gene names corresponding to the rows and columns of the regularized Laplacian

The RKHS embedding is saved to a Python dictionary with the following keys and values:
* `X` :  A NumPy array containing the values of the RKHS embedding
* `nodes`: List of nodes/gene names corresponding to the rows and columns of the RKHS embedding

**Compute HANDL embedding of target PPI network with:** `handl_embed.py`

This script computes and saves HANDL homology scores between a source and target PPI network and the HANDL embedding of the target PPI network, given the RKHS embedding of a source PPI network, the regularized Laplacian (graph kernel) of the target Network, and the list of orthologs (homologs) between the networks.

Note: other values such as the indices corresponding to Landmarks and the Landmarks used are also saved. Also, this script assumes that the given RKHS embeddings and regularized Laplacian is saved in the format outputted by `factorized_laplacian.py`.


Required parameters:

*   `-s`, `--source_rkhs_file` :  Path to source PPI network RKHS embedding
*   `-t`, `--target_laplacian_file` : Path to target PPI network regularized Laplacian
*   `-hf`, `--homologs_file` : List of homologs between source and target
*   `-o`, `--output_file` : Path to where HANDL embeddings and HANDL homology scores should be saved
*   `-n`, `--n_landmarks` : Number of landmarks HANDL should use in the embedding

Output:
The HANDL embedding is saved to a Python dictionary with the following keys and values:
*   `D` : A NumPy array containing the HANDL homology scores
*   `target_handl_C`: A NumPy array containing the values for the HANDL embedding of the target PPI network
*   `source_C`: A NumPy array containing the values for the RKHS embedding of the source PPI network
*   `target_nodes` : List of node/gene names corresponding to the rows and columns with respect to the target network and the saved matrices of embeddings/scores
*   `source_nodes` : List of node/gene names corresponding to the rows and columns with respect to the target network and the saved matrices of embeddings/scores
*   `target_landmarks` : List of node/gene names from the target used as landmarks 
*   `source_landmarks` : List of node/gene names from the source used as landmarks
*   `target_landmark_indices` : Indices of the rows of the HANDL embedding that corresponds to landmarks
*   `source_landmark_indices` : Indices of the rows of the source RKHS embedding that correspond to landmarks
*   `target_non_landmark_indices` : Indices of the rows of the HANDL embedding that correspond to homologs not used as landmarks
*   `source_non_landmark_indices` : Indices of the rows of the source RKHS embedding that correspond to homologs not used as landmarks

#### File formats
**Inputs for the scripts above:**

*   A **PPI network edge list** should be a 2 column tab separated file where each row is corresponds to an edge in the PPI network. For example, an edgelist might have a row that reads: `GENE_A GENE_B`
*   A **Homolog list** should be a 2 column tab separated file where each row corresponds to a pair of homologs between a source and a target species. For example, a Homolog list might have a row that reads: `SOURCE_GENE_A TARGET_GENE_D`

**Outputs for the scripts above:**
The Python dictionaries saved by the scripts above are saved/serialized using SciKit-Learn's JobLib module.



#### Examples
An example usage of HANDL can be found in `example/HANDL-homolog-scores` where the HANDL homology scores between proteins in fission (Sp) and baker's (Sc) yeast is computed.
You can run the example with `snakemake`, where:

*	The matrix of HANDL homology scores and HANDL embeddings with _Sc_ as the source and _Sp_ as the target will be computed and saved to `example/HANDL-homolog-scores/output/sp-to-sc-scores-and-matrices.pkl`
*	The matrix of HANDL homology scores and HANDL embeddings with _Sp_ as the source and _Sc_ as the target will be computed and saved to `example/HANDL-homolog-scores/output/sc-to-sp-scores-and-matrices.pkl`
*	The RKHS and regularized Laplacians used for HANDL embeddings are saved to `example/HANDL-homolog-scores/output/feats`.

-----
### Reference
Mark D.M. Leiserson, Jason Fan, Anthony Cannistra, Inbar Fried, Tim Lim, Thomas Schaffner, Mark Crovella, and Benjamin Hescott. (2017) "A Multi-Species Functional Embedding Integrating Sequence and Network Structure" (in submission)
