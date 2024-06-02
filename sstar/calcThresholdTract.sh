#!/bin/bash

: '
Description: Bash script to generate s* tracts for sstar scores and quantiles.
Dependencies:
    - sstar
Usage: sh simulator-Distribution-sw.sh seed num prefix_len dir_name tresh quantile_file
    - seed: seed for msprime simulation, model number
    - num: number of instances, labels with convention [0,num-1]
    - prefix_len: lenth of trees in kb
    - dir_name: directory output. MUST have following structure
        directory/
        ├─ trees/
        ├─ trees_vcf/
        ├─ egrm/
        ├─ segments/
        ├─ tar_filelist/
        ├─ sstar_score
    - thresh: float value for cut off treshold
    - quantile_file: file for quantiles summary generated with sstar quantile
Example usage:  sh simulator-Distribution-sw.sh 0 2 growth 200 ~/arg-ml/data/distribution_growth 0.5 quantile.summary.txt
Output: creates .tracts from tree files in trees_vcf and saves to directory $dir_name/$sstar_score/$seed/$thresh
Note: This step requires a QUANTILE file to be generated
Note2: the threshold results are stored in temporary directory /dev/shm which may need to be 
changed depending on the file directory set up
'

seed=$1
num=$2
prefix_len=$3
dir_name=$4
thresh=$5
quantile_file=$6

mkdir sstar_score/$seed
mkdir sstar_score/$seed/$thresh
for VAR in $(eval echo {0..$num})
do

# calculating treshold
sstar threshold --score ${dir_name}/sstar_score/200kb_${seed}_${VAR}.results  --sim-data ${quantile_file} --quantile $thresh --output /dev/shm/200kb_${seed}_${seed}_${VAR}.threshold.results

# getting tracts
sstar tract --threshold /dev/shm/200kb_${seed}_${seed}_${VAR}.threshold.results --output-prefix sstar_score/$seed/$thresh/200kb_${seed}_${SLURM_ARRAY_TASK_ID}_${VAR}.tracts

rm  /dev/shm/200kb_${seed}_${seed}_${VAR}/*

done


