import pandas as pd
# DESCRIPTION: helper script to speed up sstar accuracy calculation
for seed in tqdm([4,5,6,7,8,9,10,11,12,13,100,101,102,103,104,105,106,107,108,109]): # iterate seeds
    for model_type in ["growth", "noGrowth"]: # iterate model types
        for run in range(100): # iterate number of instances
            # load segment file
            df = pd.read_csv(f'~/arg-ml/data/distribution_{model_type}/segments/{seed}/200kb_{seed}_{run}.csv')
            # iterate number of haplotupes
            for i in range(56):
                hap1_id = 2*i
                hap2_id = 2*i+1

                segments = []
                temp = df[(df["node_sample"] == hap1_id) |(df["node_sample"] == hap2_id)][['left', 'right']]
                temp = temp.drop_duplicates()
                if temp.empty:
                    segments.append([0,0])
                else:
                    for ind in temp.index:
                        segments.append([temp['left'][ind], temp['right'][ind]])
                # create bed file for each true segment
                with open(f'~/arg-ml/data/distribution_{model_type}/segments_indv/{seed}/200kb_{seed}_{run}_{i}.bed', 'w') as o:
                    for seg in segments:
                        o.write(f'1\t{int(seg[0])}\t{int(seg[1])}\n')