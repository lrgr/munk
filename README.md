# Implementation of HANDL (Homology Assessment across Networks using Diffusion and Landmarks)

---

### From: "A Multi-Species Functional Embedding Integrating Sequence and Network Structure"

Mark D.M. Leiserson (1), Jason Fan (1), Anthony Cannistra (2), Inbar Fried (3), Tim Lim (4), Thomas Schaffner (5), Mark Crovella (4), and Benjamin Hescott (6)

* 1 Department of Computer Science and Center for Bioinformatics and Computational Biology, University of Maryland, College Park
* 2 Department of Biology, University of Washington
* 3 University of North Carolina Medical School
* 4 Department of Computer Science, Boston University
* 5 Department of Computer Science, Princeton University
* 6 College of Computer and Information Science, Northeastern University

---

### HANDL
In this repo is an implementation of HANDL. You can see an example usage in
`examples/HANDL-homolog-scores` where the HANDL homology scores between proteins in fission (Sp) and baker's (Sc) yeast is computed.
 
 You can run the example with `make all`, where the matrix of HANDL homology scores with Sc and Sp as target and source matrices will be computed and saved to disk

---

### DATA
 To run examples with our data:

##### 1 - Download and process UniProt accession ID mappings
Run `make all` in `data/name_mapping`.

##### 2 - Download and process Sc and Sp homologs from Homlogene
Run `make all` in `data/homologs`.

##### 3 - Download and process BioGRID PPI networks
Run `make all` in `data/ppi/biogrid`.
