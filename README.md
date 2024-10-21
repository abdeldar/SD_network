# SD network tools

All scripts and notebooks needed to analyze segmental duplications (SDs) from a complex networks point of view are here.

Ipython notebooks needed to reproduce results from [this paper](https://doi.org/10.1186/s12864-021-07789-7) are in **SD_network_analysis** and **Grow_network_simulations** folders.

**Grow_network_simulations**: The SD network construction, parameters fitting, analysis of overalps between SDs and CNVs and analysis of the SD network characteristics (***main_SDs_4git.ipynb***).
**Grow_network_simulations**: Simulations of network growth copying models that include UCM and PCM models and parameters fitting (***net_simul_4git.ipynb***).


All notebooks needed to reproduce results from [another publication](https://doi.org/10.1101/2023.03.18.533287) are in the following folders: ***SD_jump_models*** and ***Analyze_SD_breakpoints***.

***SD_jump_models***: The main notebook where one can build the SD network, annotate its nodes with various genomic features, run the MST-based algorithm
which predicts the minimum spanning tree of primary duplications from all edges of the SD network. All additional tests of the primary duplications prediction algorithm from [link](https://doi.org/10.1101/2023.03.18.533287) are in the ***main_second_4git.ipynb*** notebook.

NOTE: Not all input files needed for duplicated regions annotation are uploaded to the git repo. because of their size. One can either download them from the UCSC genome browser ([UCSC link]()) or skip some annotation steps or skip the whole annotation procedure and use already annotated files from the ***outputs
/annotated/*** folder.

***SDs_repeats_flanks_4git.ipynb***: Analysis of duplicated regions breakpoints, high-copy repeats and other genomic features distribution.

***NMF_Analysis***: Prediction of the duplication process signatures and signature stability tests.


## The command line tool

***predict_MST.jl***: This script allows one to construct the network of segmental duplications and distinguish primary from secondary duplications in it. The script can be run in 3 possible modes:
- without the filtering (default) (***-f 0 -k 0***).
- when all edges with matching breakpoints ("suspicious" edges) are excluded from the MST (***-f 1 -k 0***).
- when edges with matching breakpoints are excluded except for those proximal to assembly gaps (***-f 1 -k 1 --gaps_file gaps_annotation.bed***).

**To Run:**

One can run the script only if Julia is installed.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

To activate the Julia environment.
Exampl run:
```
julia --project=. predict_MST.jl -k 1 --f 1 --N=10 --input_SDs ./inputs/GRCh38GenomicSuperDup_sort.tab -g ./inputs/hg38_assembly_all_gaps_wsex.bed
```

