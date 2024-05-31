#!/bin/bash
#SBATCH --partition=chiang
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=12:00:00
#SBATCH --array=[4-13] 


eval "$(conda shell.bash hook)"
conda activate sstar

# sstar quantile --model HumanNeanderthal_4G21_constant_no_introgression.yaml --ms-dir ./msdir/ --N0 10000 --nsamp 22 --nreps 20000 --ref-index 4 --ref-size 20 --tgt-index 3 --tgt-size 2 --mut-rate 1e-8 --rec-rate 1.29e-8 --seq-len 50000 --snp-num-range 100 300 20 --output-dir quantiles_4G21  --thread 16

# sstar quantile --model HumanNeanderthal_4G21_constant_no_introgression.yaml --ms-dir ./msdir/ --N0 10000 --nsamp 22 --nreps 20000 --ref-index 4 --ref-size 20 --tgt-index 3 --tgt-size 2 --mut-rate 1e-8 --rec-rate 1.29e-8 --seq-len 50000 --snp-num-range 100 300 20 --output-dir quantiles_4G21  --thread 16
# sstar quantile --model ooa_nea_noep.yaml --ms-dir ./msdir/ --N0 10000 --nsamp 22 --nreps 10000 --ref-index 1 --ref-size 20 --tgt-index 3 --tgt-size 2 --mut-rate 2e-8 --rec-rate 1.25e-8 --seq-len 200000 --snp-num-range 400 850 25 --output-dir quantiles_SPS_200k --thread 16

# sstar quantile --model HumanNeanderthal_4G21_constant_no_introgression.yaml --ms-dir ./msdir/ --N0 10000 --nsamp 22 --nreps 1000 --ref-index 4 --ref-size 20 --tgt-index 3 --tgt-size 2 --mut-rate 1e-8 --rec-rate 1.29e-8  --seq-len 200000 --snp-num-range 550 950 25 --output-dir quantiles_4G21_200k --thread 16

# sstar quantile --model ooa_nea_noep.yaml --ms-dir ./msdir/ --N0 10000 --nsamp 224 --nreps 2000 --ref-index 1 --ref-size 112 --tgt-index 3 --tgt-size 112 --mut-rate 2e-8 --rec-rate 1.25e-8 --seq-len 50000 --snp-num-range 173 336 15 --output-dir quantiles_base --thread 16

# sstar quantile --model model_80_growth-nopulse --ms-dir ./msdir/ --N0 10000 --nsamp 56 --nreps 1000 --ref-index 2 --ref-size 28 --tgt-index 4 --tgt-size 28 --mut-rate 1e-8 --rec-rate 1.25e-8  --seq-len 200000 --snp-num-range 142 349 25 --output-dir quantiles_model_80_growth_nopulse --thread 16

# sstar quantile --model model_100_growth-nopulse --ms-dir ./msdir/ --N0 10000 --nsamp 56 --nreps 1000 --ref-index 2 --ref-size 28 --tgt-index 4 --tgt-size 28 --mut-rate 1e-8 --rec-rate 1.25e-8  --seq-len 200000 --snp-num-range 62 228 25 --output-dir quantiles_model_100_growth_nopulse --thread 16

# sstar quantile --model model_80_noGrowth-nopulse --ms-dir ./msdir/ --N0 10000 --nsamp 56 --nreps 1000 --ref-index 2 --ref-size 28 --tgt-index 4 --tgt-size 28 --mut-rate 1e-8 --rec-rate 1.25e-8  --seq-len 200000 --snp-num-range 192 380 25 --output-dir quantiles_model_80_noGrowth_nopulse --thread 16

# sstar quantile --model model_100_noGrowth-nopulse --ms-dir ./msdir/ --N0 10000 --nsamp 56 --nreps 1000 --ref-index 2 --ref-size 28 --tgt-index 4 --tgt-size 28 --mut-rate 1e-8 --rec-rate 1.25e-8  --seq-len 200000 --snp-num-range 192 380 25 --output-dir quantiles_model_100_noGrowth_nopulse --thread 16

model_type="noGrowth"
sstar quantile --model ${model_type}/model_${SLURM_ARRAY_TASK_ID}.yaml --ms-dir ./msdir/ --N0 10000 --nsamp 56 --nreps 1000 --ref-index 2 --ref-size 28 --tgt-index 4 --tgt-size 28 --mut-rate 1e-8 --rec-rate 1.25e-8  --seq-len 200000 --snp-num-range 25 450 25 --output-dir quantiles_model_${SLURM_ARRAY_TASK_ID}_${model_type} --thread 16

conda deactivate
``