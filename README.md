# Non-canonical DNA in bird T2T genomes
Code to run the analysis in the above paper (Smeds et al., 2025, bioRxiv: [https://doi.org/10.1101/2025.10.17.683159](https://doi.org/10.1101/2025.10.17.683159)), written by Linnéa Smeds.

The zebra finch T2T genome and its annotations, which are described in Formenti et al., 2025, bioRxiv ([https://doi.org/10.1101/2025.10.14.682431](https://doi.org/10.1101/2025.10.14.682431)), can be found [here](https://genomeark.s3.amazonaws.com/index.html?prefix=species/Taeniopygia_guttata/bTaeGut7/). The non-B DNA motifs annotations described in this study are uploaded [here](https://genomeark.s3.amazonaws.com/index.html?prefix=species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/nonB/).

----

Bash code is found in the directory *bash/*. Different parts of the project are separated into different bash scripts, witch descriptive names that are numbered after the order they were run. The first file ([00_setup_commands.sh](bash/00_setup_commands.sh)) contains a list of all software needed. 

Python scripts used in parts of the bash code above are found in the directory *python/*.

Code for generating the circos plots for zebra finch and chicken is found in the directory *circos/*.

All other figures for the paper are made using R, with scripts found in the directory *R/*. 

Lists of chromosome types, lengths etc can be found in the directory *helpfiles/*.
