#!/bin/bash
#SBATCH --partition=chiang
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --time=1:00:00
#SBATCH --array=[0-99] 

eval "$(conda shell.bash hook)"
conda activate sstar
level="high"
let start=${SLURM_ARRAY_TASK_ID}*10
let end=$start+10
for VAR in $(seq $start $end)
do
sstar threshold --score data/detectSite/trees_${level}/sstar_scores_50kb/10kb_42_${VAR}.results --sim-data sstar/quantiles_test_1k/quantile.summary.txt  --quantile 0.99 --output data/detectSite/trees_${level}/sstar_thresholds_50kb/50kb_42_${VAR}.threshold.results
done
conda deactivate
