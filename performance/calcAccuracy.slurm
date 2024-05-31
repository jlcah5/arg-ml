#!/bin/bash
#SBATCH --partition=chiang
#SBATCH --nodes=1 
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:00:00
#SBATCH --array=[12] 


# python calcAccuracy_nn.py Distribution transformer noGrowth ${SLURM_ARRAY_TASK_ID} 50000
# python calcAccuracy_nn.py Distribution transformer noGrowth ${SLURM_ARRAY_TASK_ID} 250000
# python calcAccuracy_nn.py Distribution transformer noGrowth ${SLURM_ARRAY_TASK_ID} 600000

# python calcAccuracy_nn.py Distribution transformer growth ${SLURM_ARRAY_TASK_ID} 50000
# python calcAccuracy_nn.py Distribution transformer growth ${SLURM_ARRAY_TASK_ID} 250000
python calcAccuracy_nn.py Distribution transformer growth ${SLURM_ARRAY_TASK_ID} 600000


# mkdir /scratch1/jcahoon/argml/data/distribution_growth_validate/segments_indv/${SLURM_ARRAY_TASK_ID}
# mkdir /scratch1/jcahoon/argml/data/distribution_noGrowth_validate/segments_indv/${SLURM_ARRAY_TASK_ID}
# python sstar_segindv.py ${SLURM_ARRAY_TASK_ID}
# python calcAccuracy_sstar.py growth ${SLURM_ARRAY_TASK_ID}
# python calcAccuracy_sstar.py noGrowth ${SLURM_ARRAY_TASK_ID}

# python calcAccuracy_nn.py SPS transformer transfer ${SLURM_ARRAY_TASK_ID}
# python calcAccuracy_nn.py Gower transformer transfer ${SLURM_ARRAY_TASK_ID}

#python calcAccuracy_nn.py SPS transformer infer ${SLURM_ARRAY_TASK_ID}
# python calcAccuracy_nn.py SPS transformer true ${SLURM_ARRAY_TASK_ID} # switch cut off
# python calcAccuracy_nn.py SPS transformer eval ${SLURM_ARRAY_TASK_ID}
# python calcAccuracy_nn.py Gower transformer infer ${SLURM_ARRAY_TASK_ID}
# python calcAccuracy_nn.py Gower transformer true  ${SLURM_ARRAY_TASK_ID} # switch cut off
# python calcAccuracy_nn.py Gower transformer eval ${SLURM_ARRAY_TASK_ID}


# python calcAccuracy_nn.py SPS transformer none
# python calcAccuracy_nn.py SPS cnn none
# python calcAccuracy_nn.py Gower transformer none
# python calcAccuracy_nn.py Gower transformer 08
# python calcAccuracy_nn.py Gower transformer Bal

# python calcAccuracy_nn.py Gower cnn none
