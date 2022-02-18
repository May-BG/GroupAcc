# set -uex

#!/bin/bash

#PBS -l nodes=1:ppn=1

#PBS -l walltime=10:00:00

#PBS -l mem=20G

#PBS -A yuh371_a_g_gc_default

# Get started

# Go to the correct place
module use /gpfs/group/RISE/sw7/modules

module load r/4.0.3

cd /storage/home/xmz5176/group/data/hars/2013tfbs/processed/final/filter1_50_zip/bed/noover/filtered
f=$1
cd $f
#cd /storage/home/xmz5176/group/data/hars/2013tfbs/phangorn_110/single_nogap/polr3

Rscript /storage/home/xmz5176/group/data/hars/code/mixSingle.R $f /storage/home/xmz5176/group/data/hars/2013tfbs/tfbs_filtered/merged_tfbs2013/info5/tfbs_ref_pml.RData /storage/home/xmz5176/group/data/hars/2013tfbs/processed/final/tfbs_n1_nogap.nh hg19,rheMac2 ${f}/pip > ${f}/pip/mixSingle.txt

