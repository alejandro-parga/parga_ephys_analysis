import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

#hard coded variables for testing
# file = "analysis/NHP/VITable_merged.csv"
# outpath = "analysis/figures/NHP/VITable/"
# outfile = "analysis/figures/NHP/VITable/done.txt"

file = snakemake.input.merged
outfile = snakemake.output.out
outpath = os.path.dirname(outfile)

#check if output directory exists
#if not, create
isExist = os.path.exists(outpath)
if not isExist:
    print(f"creating output directory: {outpath}")
    os.makedirs(outpath)

data = pd.read_csv(file)

colnames = list(data.columns)

cells = data['cell_name']

for c in colnames[1:]:
    plt.figure(constrained_layout=True)
    plt.scatter(cells,data[c])
    plt.xlabel('cell name')
    plt.ylabel(c)
    plt.xticks(rotation=45,ha="right")
    plt.title(c)
    plt.axhline(y=np.nanmean(data[c]),color="black",linestyle="--",label="mean")
    plt.savefig(f'{outpath}/{c}.png')
    plt.close()

f=open(outfile, "w")
f.write("Done!")