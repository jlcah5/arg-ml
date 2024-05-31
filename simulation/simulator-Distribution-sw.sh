#!/bin/bash
: '
Description: Bash script to generate validation data
Dependencies:
    - python3: pandas, tskit, demes, tqdm
    - relate
    - egrm
Usage: sh simulator-Distribution-sw.sh seed num length model_type directory window_size step
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
        ├─ tar_filelist/
    - window_size: size of window in base pairs
    - step: size of step size in base pairs
Example usage:  sh simulator-Distribution-sw.sh 0 2 200000 growth /scratch1/jcahoon/arg-ml/data/distribution_growth/ 50000 10000
Output: creates subdirectories in each {trees, trees_vcf, egrm, segments} labeled $seed containing ($num-$window_size)/$step instances for each file type
    - trees: tree file
    - trees_vcf: vcf file generated from tree file
    - egrm: npy files for egrm matric for both true and inferred trees (size: [112,112]) and label (size: [,112])
    - segments: true segments in csv files
    - tar_filstlist: tarballs of neural network inputs
'

seed=$1
let num=$2
let length=$3
model_type=$4
dir_name=$5
let window_size=$6
let step=$7

let window_number=$length-${window_size}
let window_number=${window_number}/$step
let length_kb=$length/1000

#### Data Generation ####
python simulator-Distribution-sw.py $seed $num ${length_kb} $seed ${model_type} ${dir_name} ${window_size} ${step}
########################

#### ARG inference ####
change to temporary directory for file I/O
# may need to change to appropriate directory
cd /dev/shm/

let num=$num-1 # minus one since range is inclusive [0,num]
for i in $(eval echo {0..$num})
do
for j in $(eval echo {0..${window_number}})
do
let window=${j}*10000
prefix="${length_kb}kb_${seed}_${i}_${window}"

# convert from vcf
RelateFileFormats --mode ConvertFromVcf --haps ${prefix}.haps --sample ${prefix}.sample -i ${dir_name}/trees_vcf/${seed}/${length_kb}kb_${seed}_${i}_${window}

# remove non biallelic snps
RelateFileFormats --mode RemoveNonBiallelicSNPs --haps ${prefix}.haps -o ${prefix}.clean

# infer
Relate --mode All -m 1.25e-8 -N 20000 --haps ${prefix}.clean.haps --sample ${prefix}.sample --map ~/genetic_map.txt --seed 1 -o ${prefix}

# convert to tree sequence
RelateFileFormats --mode ConvertToTreeSequence -i ${prefix} -o ${prefix}.infer

# convert to egrm
trees2egrm --output-format numpy ${prefix}.infer.trees  --c --haploid --output ${length_kb}kb_${seed}_${i}_${window}.infer.input

# clean files
mv ${length_kb}kb_${seed}_${i}_${window}.infer.input.npy ${dir_name}/egrm/${seed}/${length_kb}kb_${seed}_${i}_${window}.infer.input.npy
rm ${length_kb}kb_${seed}_${i}_${window}.infer.input_mu.npy ${length_kb}kb_${seed}_${i}_${window}.infer.input.log.log
rm -rf ${prefix}.*
done
done
# ########################


#### Tarball creation ####

# clearing files
rm ${dir_name}/tar_filelist/file_${seed}.txt
touch ${dir_name}/tar_filelist/file_${seed}.txt

# creating file list
for i in $(eval echo {0..$num})
do
for j in $(eval echo {0..${window_number}})
do
let window=${j}*10000
echo ${dir_name}/egrm/${seed}/${length_kb}kb_${seed}_${i}_${window}.input.npy >> ${dir_name}/tar_filelist/file_${seed}.txt
echo ${dir_name}egrm/${seed}//${length_kb}kb_${seed}_${i}_${window}.infer.input.npy >> ${dir_name}/tar_filelist/file_${seed}.txt
echo ${dir_name}/egrm/${seed}/${length_kb}kb_${seed}_${i}_${window}.output.npy >> ${dir_name}/tar_filelist/file_${seed}.txt
done
done

# zipping files together
tar -c -T ${dir_name}//tar_filelist/file_${seed}.txt -f ${dir_name}/tar_filelist/200kb_${seed}.tar
echo "finished"
