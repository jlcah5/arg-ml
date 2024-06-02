#!/bin/bash
: '
Description: Bash script to generate s* scores
Dependencies:
    - sstar
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
Note: If changing simulation from 56 haplotypes, please change ref.ind.lisr and tgt.ind.list to 
list the haplotypes of the correct reference and target haplotypes, respectively.
Note 2: sstar defaults step and window size as 10kb and 50kb, respectively. If changing these 
parameters, make sure to update how the sstar score is calculated.
'
let seed=$1
num=$2
prefix_len=$3
dir_name=$4
mkdir ${dir_name}/sstar_score
for VAR in $(eval echo {0..$num})
do
    sstar score --vcf ~ ${dir_name}/sstar_score/trees_vcf/${seed}/${prefix_len}kb_${seed}_${VAR}.vcf  --ref ref.ind.list --tgt tgt.ind.list --output   ${dir_name}/sstar_score/${prefix_len}kb_${seed}_${VAR}.results 
done