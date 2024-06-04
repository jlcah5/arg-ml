# Directory Contents
## ML
Contents: model class files, trained models, training scripts, model prediction
Important scripts:
- indv-50kb-distribution-infer.py: trains model and saves weights
- indv-50kb-distribution-sw-infer.py: evaluates trained model on sliding windows

## Performance
Contents: scripts for accuracy calculation
Important scripts:
- calcAccuracy_nn.py: calculates files necessary to generate precision-recall curves for ml benchmark
- calcAccuracy_sstar.py: calculates files necessary to generate precision-recall curves for sstar benchmark

## Plot
Contents: notebooks to generate precision recall curves, chromosome paintings, and demographic model plots
Important scripts:
- distributionDemographies.ipynb: generate demographic model plots
- precision-recall.ipynb: generate pr curves and chromosome paintings
## Simulation
Contents: scripts to simulate data, sampled demographic models 
- simulator-Distribution.sh: generate training data
- simulator-Distribution-sw.sh: generate evaluation data with sliding window
## Sstar
Contents: scripts to run sstar benchmark
- calcSstar.sh: calculates s* scores
- quantile.sh: generates null distribution data
- calcThresholdTract.sh: identifies introgressed segments
## Notes
- If using the scripts directory, make sure to update the paths
- Some data files are too big to upload over github and will need to be transferred another way.
- Data generation requires a specific file structure. Please organize the directory in this structure:
	-	directory/
        ├─ trees/
        ├─ trees_vcf/
        ├─ egrm/
        ├─ segments/
        ├─ tar_filelist/
        ├─ segments_indv/
- The calcThresholdTract.sh can only be run AFTER calcSstar.sh and quantile.sh are completed.
- If the simulation parameters are altered (like length of tree series), other scripts may need to be updated
