#!/bin/bash
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=1:00:00
#SBATCH --array=[4-13] 


# sstar score --vcf test.data.vcf --ref test.ref.ind.list --tgt test.tgt.ind.list --output test.score.results

# let start=${SLURM_ARRAY_TASK_ID}*100
# let end=$start+99
# for VAR in $(seq $start $end)

let seed=${SLURM_ARRAY_TASK_ID}
model="noGrowth"
mkdir /scratch1/jcahoon/argml/data/distribution_${model}_validate/sstar_score
for VAR in $(seq 0 99)
do
#     sstar score --vcf  data/detectSite/trees_SPS_validateSstar/trees_vcf/50kb_0_${VAR}.vcf --ref ref.SPS.ind.list --tgt tgt.SPS.ind.list --output data/detectSite/trees_SPS_validateSstar/sstar_scores/50kb_0_${VAR}.results
#     sstar score --vcf  data/detectSite/trees_Gower_validateSstar/trees_vcf/200kb_42_${VAR}.vcf --ref ref.ind.list --tgt tgt.ind.list --output data/detectSite/trees_Gower_validateSstar/sstar_scores/200kb_42_${VAR}.results 
#     sstar score --vcf  data/detectSite/trees_SPS_validateSstar/trees_vcf/200kb_42_${VAR}.vcf --ref ref.ind.list --tgt tgt.ind.list --output data/detectSite/trees_SPS_validateSstar/sstar_scores/200kb_42_${VAR}.results 
    sstar score --vcf /scratch1/jcahoon/argml/data/distribution_${model}_validate/trees_vcf/${seed}/200kb_${seed}_${VAR}.vcf  --ref ref.ind.list --tgt tgt.ind.list --output  /scratch1/jcahoon/argml/data/distribution_${model}_validate/sstar_score/200kb_${seed}_${VAR}.results 
done