import pandas as pd
import numpy as np
import warnings
import time
import os

feat = pd.read_csv("weighted_feat2.csv")
net = pd.read_csv(os.path.join("clean_data", "net.csv"))

print("Files read. Time for computation!", flush=True)

pd.options.mode.chained_assignment = None
t0 = time.time()
i = 0

for d in net.itertuples():
    if i != d.chr_num:
        i = d.chr_num
        print(f"Now on chromosome: {i}", flush=True)
        t = time.time() - t0
        print(f"Time is {t//3600} hours {(t%3600) // 60} minutes {(t%3600) % 60} seconds", flush=True)

    condition1 = (feat["strand"] == d.strand)
    condition2 = (feat["chr_num"] == d.chr_num) 
    condition3 = (feat["left"] < d.pos)
    condition4 = (feat["right"] > d.pos)
    df = feat[condition1 & condition2 & condition3 & condition4]
    if len(df) != 0:
        net.at[d.Index, "total"] = df["reads"].sum()
        net.at[d.Index, "feature"] = df[df["reads"] == max(df["reads"])]["class"].iloc[0]
        
        
print("Job done. Now making the file", flush=True)
net.to_csv("weighted_net.csv")