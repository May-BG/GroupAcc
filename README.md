# GroupAcc

## Description
The GroupAcc includes two pooling-based phylogenetic approaches with improved statistical power to infer weakly accelerated evolution. The key idea of GroupAcc is to group TFBSs by the bound transcription factor and then examine whether each TFBS group as a whole shows an elevated substitution rate in a lineage of interest.

## Requirements
The GroupAcc is implemented in R with several packages.

### Input files
X.zip: a zip directory composed of alignments of specific group of genetic elements

tfbs_ref_pml.RData: a reference phylogenetic model file in the format of RData

tfbs_n1_nogap.nh: a reference phylogenetic tree containing the topology

hg19: the lineage of interest

### GroupLRT
``` Rscript groupLRT_1221.R ${h}.zip tfbs_ref_pml.RData tfbs_n1_nogap.nh hg19 $f > ${f}/hg19.txt```
```
mkdir groupLRT_output/
Rscript groupLRT.R ZZZ3.zip tfbs_ref_pml.RData tfbs_n1_nogap.nh hg19 groupLRT_output/ > groupLRT_output/hg19_groupLRT.txt
```
### Mixture Model
```$ Rscript mixSingle_1221.R BDP1 tfbs_ref_pml.RData tfbs_n1_nogap.nh hg19,panTro2 ${f}/panTro2 > ${f}/panTro2/mixSingle.txt```
```
mkdir mixSingle_output/
Rscript mixSingle.R POU5F1/ tfbs_ref_pml.RData tfbs_n1_nogap.nh hg19,panTro2 mixSingle_output/ > mixSingle_output/mixSingle.txt
```

## Help

## Authors
