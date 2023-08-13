import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import seaborn as sns

#hard coded variables for testing
# animal = "Etv1"
# cell = "AP030923B"
# outfile = f"analysis/figures/{animal}/FITable/done.txt"
# outpath = os.path.dirname(outfile)

animal = snakemake.params.animal
outfile = snakemake.output.out
outpath = os.path.dirname(outfile)

#check if output directory exists
#if not, create
isExist = os.path.exists(outpath)
if not isExist:
    print(f"creating output directory: {outpath}")
    os.makedirs(outpath)

print(f"animal: {animal}")
print(f"outfile = {outfile}")
print(f"outpath = {outpath}")

#get file name base on animal
file = f"analysis/{animal}/FITable_merged.csv"

#read file
data = pd.read_csv(file)

#get cellnames
cellnames = data["cell_name"]
#get data and plot for each cell
for cell in cellnames:
    #initialize empty data frame
    df = pd.DataFrame()

    print(cell)
    #subsetdata for cell
    subdata = data[data.cell_name == cell]
    subdata.reset_index(drop=True,inplace=True)
    #get AP data
    APstr=subdata["APs[][0]"].astype(str)
    #remove brackets
    APstr1=APstr[0].replace("[","").replace("]","")
    #split numbers into list
    APstr2=re.split(" ",APstr1)
    #convert character into integer
    APnum=[int(n) for n in APstr2]
    #add values to dataframe
    df["AP"]=APnum
    #get current injection
    CIstr=subdata["current_injection[][0]"].astype(str)
    CIstr1=CIstr[0].replace("[","").replace("]","")
    CIstr2=re.split(" ",CIstr1)
    CInum=[int(n) for n in CIstr2]
    df["CI"]=CInum

    #plot current injection vs APs
    plt.figure(constrained_layout=True)
    plt.scatter(x=df["CI"],y=df["AP"])
    plt.title(cell)
    plt.xlabel("Current Injection")
    plt.ylabel("Number of APs")
    plt.savefig(f'{outpath}/{cell}.png')
    plt.close()

f=open(outfile, "w")
f.write("Done!")
