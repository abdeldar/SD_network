# SD network tools

All scripts and notebooks needed to analyze segmental duplications (SDs) from a complex network point of view are here.

## The command line tool ##

***predict_MST.jl***:      This script allows one to construct the network of segmental duplications and distinguish primary from secondary duplications in it. The script can be run in 3 possible modes: <br />
- without the filtering (default) (***-f 0 -k 0***).
- when all edges with matching breakpoints ("suspicious" edges) are excluded from the MST (***-f 1 -k 0***).
- when edges with matching breakpoints are excluded except for those proximal to assembly gaps (***-f 1 -k 1 --gaps_file gaps_annotation.bed***).

**Quick Run:**

One can run the script only if Julia is installed.

To activate the Julia environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Example runs:

```
julia --project=. predict_MST.jl --help

julia --project=. predict_MST.jl -k 0 --f 0 --input_SDs ./inputs/GRCh38GenomicSuperDup_sort.tab

julia --project=. predict_MST.jl -k 1 --f 1 --N=10 --input_SDs ./inputs/GRCh38GenomicSuperDup_sort.tab -g ./inputs/hg38_assembly_all_gaps_wsex.bed
```

Only ***-i*** argument which points to the list of SDs file is obligatory. One can find the list of _hg38_ segmental duplications at ***inputs/GRCh38GenomicSuperDup_sort.tab***.

NOTE: The assembly gap annotation for the _hg38_ can be found at ***inputs/hg38_assembly_all_gaps_wsex.bed***. For other gemnomes use another input annotations.

 
 
 
## Jupyter notebooks ##

Jupyter notebooks needed to reproduce results from [this paper](https://doi.org/10.1186/s12864-021-07789-7) are in **SD_network_analysis** and **Grow_network_simulations** folders.

**Grow_network_simulations**: The SD network construction, parameters fitting, analysis of overalps between SDs and CNVs and analysis of the SD network characteristics (***main_SDs_4git.ipynb***).

**Grow_network_simulations**: Simulations of network growth copying models (UCM and PCM) and parameters fitting (***net_simul_4git.ipynb***).


All jupyter notebooks needed to reproduce results from [another publication](https://doi.org/10.1101/2023.03.18.533287) are in the following folders: ***SD_jump_models*** and ***Analyze_SD_breakpoints***.

***SD_jump_models***: Includes the main notebook where one can build the SD network, annotate its nodes with various genomic features, run the MST-based algorithm which predicts the minimum spanning tree of primary duplications from all edges of the SD network. All additional tests of the MST-based algorithm from [link](https://doi.org/10.1101/2023.03.18.533287) paper are in the ***main_second_4git.ipynb*** notebook.

NOTE: Input files needed for duplicated regions annotation are NOT uploaded to the git repository because of their size. One can either download them from the UCSC genome browser ([UCSC link](https://genome.ucsc.edu/cgi-bin/hgTables)) or skip some annotation steps or skip the whole annotation procedure and use already pre-annotated files from the ***outputs/annotated/*** folder.

***SDs_repeats_flanks_4git.ipynb***: Analysis of duplicated regions breakpoints, high-copy repeats and other genomic features distribution.

***NMF_Analysis***: Prediction of the duplication process signatures and signature stability tests.

