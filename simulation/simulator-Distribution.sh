#!/bin/bash
: '
Description: Bash script to generate training data
Dependencies:
    - python3: pandas, tskit, demes, tqdm
    - relate
    - egrm
Usage: sh simulator-Distribution.sh seed num length model_type directory
    - seed: seed for msprime
    - num: number of instances, labels with convention [0,num-1]
    - length: lenth of trees in base pairs
    - model_type: choices {growth, noGrowth}
    - dir_name: directory output. MUST have following structure
        directory/
        ├─ trees/
        ├─ trees_vcf/
        ├─ egrm/
        ├─ segments/
Example usage: sh simulator-Distribution-sw.sh 0 1 200000 growth ~/arg-ml/data/distribution_growth/ 50000 10000
Output: creates subdirectories in each {trees, trees_vcf, egrm, segments} labeled with $seed containing $num instances for each file type.
    - trees: tree file
    - trees_vcf: vcf file generated from tree file
    - egrm: npy files for egrm matric for both true and inferred trees (size: [112,112]) and label (size: [,112])
    - segments: true segments in csv files
'

# user input
seed=$1
let num=$2 
let length=$3/1000 
model_type=$4 
dir_name=$5 

python simulator-Distribution.py $seed $num $length $seed ${model_type} ${dir_name} 

# change to temporary directory for file I/O
## may need to change to appropriate directory
cd /dev/shm/

# infer arg with relate
let num=$num-1 # minus one since range is inclusive [0,num]
for i in $(eval echo {0..$num})
do
prefix="${length}kb_${seed}_${seed}_${i}"

# convert from vcf
RelateFileFormats --mode ConvertFromVcf --haps ${prefix}.haps --sample ${prefix}.sample -i ${dir_name}/trees_vcf/${seed}/${length}kb_${seed}_${seed}_${i}

# remove non biallelic snps
RelateFileFormats --mode RemoveNonBiallelicSNPs --haps ${prefix}.haps -o ${prefix}.clean

# infer
Relate --mode All -m 1.25e-8 -N 20000 --haps ${prefix}.clean.haps --sample ${prefix}.sample --map genetic_map.txt --seed 1 -o ${prefix}

# convert to tree sequence
RelateFileFormats --mode ConvertToTreeSequence -i ${prefix} -o ${prefix}.infer

# convert to egrm
trees2egrm --output-format numpy ${prefix}.infer.trees  --c --haploid --output ${length}kb_${seed}_${seed}_${i}.infer.input

# clean output files
mv ${length}kb_${seed}_${seed}_${i}.infer.input.npy ${dir_name}/egrm/${seed}/${length}kb_${seed}_${seed}_${i}.infer.input.npy
rm ${length}kb_${seed}_${seed}_${i}.infer.input_mu.npy ${length}kb_${seed}_${seed}_${i}.infer.input.log.log
rm -rf ${prefix}.*
done

echo "finished"
#clean up the vcf files
rm -rf /scratch1/jcahoon/argml/data/distribution_${model_type}/trees_vcf/${seed}
