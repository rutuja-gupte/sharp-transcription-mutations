import pandas as pd
import numpy as np
import warnings
import time
import os

feat = pd.read_csv(os.path.join("clean_data", "features.csv"))
net = pd.read_csv(os.path.join("clean_data", "net.csv"))
feat["reads"] = 0

print("Files read. Time for computation!")

pd.options.mode.chained_assignment = None
t0 = time.time()
i = 0

for d in net.itertuples():
    if i != d.chr_num:
        i = d.chr_num
        print(f"Now on chromosome: {i}")
        t = time.time() - t0
        print(f"Time is {t//3600} hours {(t%3600) // 60} minutes {(t%3600) % 60} seconds")

    condition1 = (feat["strand"] == d.strand)
    condition2 = (feat["chr_num"] == d.chr_num) 
    condition3 = (feat["left"] < d.pos)
    condition4 = (feat["right"] > d.pos)
    df = feat[condition1 & condition2 & condition3 & condition4]
    dists = list(1/(d.pos - df["left"]))
    total = sum(dists)
    for j in range(len(df)):
        feat.at[(df.iat[j,0]-1), "reads"] += d.reads * dists[j] / total

print("Job done. Now making the file")
feat.to_csv("weighted_feat.csv")