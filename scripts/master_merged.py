import csv
import os
import glob
import pandas as pd

# assign variable animal
animal = snakemake.params.animal
outpath = snakemake.output.merged


# hard coded variables for testing
# print("TESTING")
# animal = "Etv1"
# outpath = f"analysis/{animal}/membrane_properties_{animal}.csv"

# for log print variables
print(f"animal: {animal}")
print(f"output file: {outpath}")

#check if output directory exists
#if not, create
outdir = os.path.dirname(outpath)
isExist = os.path.exists(outdir)
if not isExist:
    print(f"creating output directory: {outdir}")
    os.makedirs(outdir)

# get the files for this animal
FIdata = pd.read_csv(f"analysis/{animal}/FITable_merged.csv")
VIdata = pd.read_csv(f"analysis/{animal}/VITable_merged.csv")
ZAPdata = pd.read_csv(f"analysis/{animal}/ZAPTable_merged.csv")  
APdata = pd.read_csv(f"analysis/{animal}/APTable_merged.csv")                   

# merge data
merge1 = pd.merge(FIdata,VIdata, how="outer", on="cell_name")
merge2 = pd.merge(merge1,ZAPdata, how="outer", on="cell_name")
merge3 = pd.merge(merge2,APdata, how="outer", on="cell_name")
merge3.to_csv(outpath, index=False)

