# Homology Assessment across Networks using Diffusion and Landmarks (HANDL)

### HANDL
HANDL is an algorithm for embedding proteins in a target network (e.g. mouse) into a source network (e.g. from human). The HANDL algorithm was developed by Mark Crovella (Boston University), Benjamin Hescott (Northeastern University), and Max Leiserson (University of Maryland, College Park), and their respective research groups.

### Setup
#### Dependencies
**Required Programs:** GNU Make, wget, Python 2.7.12+, PIP
**Python Packages:** Install required modules with `pip install -r requirements.txt`

#### Data
HANDL requires as input a source and target PPI network, and a list of homologs mapping a subset of the nodes in each network. We include scripts for downloading and processing data for _S. cerevisiae_ (_Sc_) and _S. pombe_ (_Sp_) networks and homologs, and mapping them into the same namespace.

You can process the _Sc_ and _Sp_ data with the following commands.

1. **Download and process UniProt accession ID mappings.** Run `make all` in `data/name_mapping`.

2. **Download and process Sc and Sp homologs from Homlogene.** Run `make all` in `data/homologs`.

3. **Download and process BioGRID PPI networks.** Run `make all` in `data/ppi/biogrid`.

### Usage

#### Scripts and command-line arguments
**Computing RKHS factors/factorization of regularized laplacian.** `factorized_laplacian.py`

**Computing HANDL embedding.** `handl_embed.py`

#### File formats
Input/Output for the above scripts.

#### Examples
An example usage of HANDL can be found in `examples/HANDL-homolog-scores` where the HANDL homology scores between proteins in fission (Sp) and baker's (Sc) yeast is computed.
You can run the example with `make all`, where the matrix of HANDL homology scores with _Sc_ and _Sp_ as target and source matrices will be computed and saved to `/path/to/output/`.

### Reference
Mark D.M. Leiserson, Jason Fan, Anthony Cannistra, Inbar Fried, Tim Lim, Thomas Schaffner, Mark Crovella, and Benjamin Hescott. (2017) "A Multi-Species Functional Embedding Integrating Sequence and Network Structure" (in submission)
