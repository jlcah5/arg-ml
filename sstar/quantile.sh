#!/bin/bash

: '
Description: Bash script to generate quantile files
Dependencies:
    - sstar
    - ms
Usage: sh simulator-Distribution-sw.sh seed num prefix_len dir_name
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
Example usage:  sh simulator-Distribution-sw.sh 0 2 growth 200 ~/arg-ml/data/distribution_growth
Output: creates .results from tree files in trees_vcf and saves to directory dir_name/sstar_score
Note: Update --ms-dir to correct path
Note2: May need to update -N0, -nsamp, --nreps, --ref-index, --ref-size, --tgt-index, --tgt-size
--mut-rate --rec-rate, --snp-num-range etc if simulation is altered
'
seed=$1
model=$2
length=$3
sstar quantile --model ~/arg-ml/simulation/demographic_models/${model}/model_${seed}.yaml --ms-dir ./msdir/ --N0 10000 --nsamp 56 --nreps 1000 --ref-index 2 --ref-size 28 --tgt-index 4 --tgt-size 28 --mut-rate 1e-8 --rec-rate 1.25e-8  --seq-len $length --snp-num-range 25 450 25 --output-dir quantiles_model_${seed}_${model} --thread 16

conda deactivate
``