#!/bin/bash
#SBATCH --partition=main
#SBATCH --nodes=1 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2:00:00

module load htslib
module load bcftools
"""
USAGE: sh gnomadTabix.sh Chr Base Max
    Chr: chromosome
    Base: starting window of the chromosome
    Max: number of windows 
See suggestions for chr 20-22 below
TO-DO: 
- Fill out PATH_TO_DATA/ PATH_TO_OUTPUT
- adjust relate simulation
- create tarball for input into model
"""
chr=$1
base=$2 
max=$3

# 20
# max=16
# base=500000

# 21
# max=9
# base=5500000

# 22
# max=9
# base=10500000

zero=0;
for i in $( eval echo {0..${max}} )
do
let start=$base+4000000*$i
let end=$start+5000000

ret=$(tabix /project/jitang_1167/data/1KG+HGDP_gnomADr3.1.2_phased/phased_haplotypes_v2_vcf/hgdp1kgp_chr$chr.filtered.SNV_INDEL.phased.shapeit5.vcf.gz chr$chr:$start-$end | wc -l);

if [ $ret -eq $zero ]
then
    echo "false";
    exit;
else
    # extract window
    tabix -h PATH_TO_DATA/hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.vcf.gz chr$chr:$start-$end | bgzip > hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.$start.$end.vcf.gz

    # for every window, extract the population
    for pop_file in ~/argml/gnomad/gnomad_subpops/*
    do 
    
    pop=${pop_file%.txt}
    pop=${pop##*/}
    prefix="${pop}_chr${chr}_${start}_${end}_${i}"
    bcftools view -S ${pop_file} hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.$start.$end.vcf.gz | bgzip > ${prefix}.vcf.gz
        
    # infer arg
    RelateFileFormats --mode ConvertFromVcf --haps ${prefix}.haps --sample ${prefix}.sample -i ${prefix}
    # remove non biallelic snps
    RelateFileFormats --mode RemoveNonBiallelicSNPs --haps  ${prefix}.haps -o ${prefix}.clean
    # infer tree
    Relate --mode All -m 1.25e-8 -N 20000 --haps ${prefix}.clean.haps --sample ${prefix}.sample --map ~/genetic_map.txt --seed 1 -o ${prefix}
    # convert to tree sequence
    RelateFileFormats --mode ConvertToTreeSequence -i ${prefix} -o ${prefix}.infer
    # calculate GRM
    trees2egrm --output-format numpy ${prefix}.infer.trees --c --haploid --output PATH_TO_OUTPUT/$prefix
    # clean  
    rm -rf ${prefix}.*

    done

fi

done
# to-do: tarball of PATH_TO_OUTPUT/*.npy files

