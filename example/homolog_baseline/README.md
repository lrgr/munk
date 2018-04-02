# SL prediction baseline with homologs

In this baseline, we say that for a pair of genes (a,b) in a source species, and a pair of genes (u,v) in a target species, that (u,v) is SL iff (a,b) is SL and (a,u) and (b,v) are homologs.

### Notes:
- Inconclusive GIs are excluded
- For BioGrid, non SLs are sampled from cartesian product of nodes obtained from PPI networks
- `-s` option for `homolog_mapping_table.py` triggers the option to sample non SLs

### Usage:
Run 'Snake --configfile configs/{collins-roguev.yml | biogrid.v3.4.157}' to obtain results for Collins/Roguev or BioGrid v3.4.157.

### Results:

| Source Species       | Target Species       | Precision | Recall   | FPR      | F1       | 
|----------------------|----------------------|-----------|----------|----------|----------| 
| Collins (Sc)         | Roguev (Sp)          | 0.269231  | 0.002493 | 0.000384 | 0.004940 | 
| Roguev (Sp)          | Collins (Sc)         | 0.291667  | 0.001934 | 0.000271 | 0.003842 | 
| Biogrid 3.4.157 (Sc) | Biogrid 3.4.157 (Sp) | 1.000000  | 0.005203 | 0.000000 | 0.010352 | 
| Biogrid 3.4.157 (Sp) | Biogrid 3.4.157 (Sc) | 1.000000  | 0.000603 | 0.000000 | 0.001204 | 

