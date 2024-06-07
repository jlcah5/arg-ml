#!/bin/bash
: '
    Description: Bash script to simulate long tree sequences, infer the arg, then window the data
Dependencies:
    - python3: pandas, tskit, demes, tqdm
    - relate
    - egrm
Usage: sh simulator-Distribution-sw.sh seed num length model_type directory window_size step
    - seed: seed for msprime
    - num: number of instances, labels with convention [0,num-1]
    - length: lenth of trees in kb
    - model_type: choices {growth, noGrowth}
    - dir_name: directory output. MUST have following structure
        directory/
        ├─ trees/
        ├─ trees_vcf/
        ├─ egrm/
        ├─ segments/
        ├─ tar_filelist/
    - window_size: size of window in base pairs
    - step: size of step size in base pairs
Example usage:  sh simulator-Distribution-sw-long.sh 0 2 2000 growth /scratch1/jcahoon/arg-ml/data/distribution_growth/ 50000 10000
Output: creates subdirectories in each {trees, trees_vcf, egrm, segments} labeled $seed containing ($num-$window_size)/$step instances for each file type
    - trees: tree file
    - trees_vcf: vcf file generated from tree file
    - egrm: npy files for egrm matric for both true and inferred trees (size: [112,112]) and label (size: [,112])
    - segments: true segments in csv files
    - tar_filstlist: tarballs of neural network inputs
TO-DO: tarball creation -- to vary overhang
'

seed=$1
let num=$2
let length=$3
model_type=$4
dir_name=$5
let window_size=$6
let step=$7

# generate the long segments
python simulator-Distribution-sw-long.py $seed $num $length $seed ${model_type} ${dir_name}

# infer the arg
for i in $(eval echo {0..$num})
do
prefix="${length}kb_${seed}_${i}"
# convert from vcf
echo "convert to vcf: ${dir_name}/trees_vcf/${base_seed}/${length}kb_${base_seed}_${i}"
RelateFileFormats --mode ConvertFromVcf --haps ${prefix}.haps --sample ${prefix}.sample -i ${dir_name}/trees_vcf/${seed}/${length}kb_${seed}_${i}

# remove non biallelic snps
RelateFileFormats --mode RemoveNonBiallelicSNPs --haps ${prefix}.haps -o ${prefix}.clean
echo "convert to infer"
# infer
Relate --mode All -m 1.25e-8 -N 20000 --haps ${prefix}.clean.haps --sample ${prefix}.sample --map ~/genetic_map.txt --seed 1 -o ${prefix}
echo "convert to trees"
# convert to tree sequence
RelateFileFormats --mode ConvertToTreeSequence -i ${prefix} -o ${dir_name}/trees/${seed}/${prefix}.infer

rm -rf ${prefix}.*

done
echo "finished"

# window the tree
python simulator-Distribution-sw-long-split.py $seed $num $length $seed ${model_type}  ${dir_name} ${window_size} ${step}


#### Tarball creation ####
