import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import seaborn as sns
import glob

# hard coded variables for testing
# print ("TESTINGS")
# animal = "NHP" #options:NHP, Thy1

animal = snakemake.params.animal

# get list of cells
# cellnames = os.listdir(f"results/{animal}/")
filepath = f"results/{animal}/*/*_0noiseF_Table.csv"
files = glob.glob(filepath)
dirs = [os.path.dirname(f) for f in files]
cellnames = [d.replace(f"results/{animal}/","") for d in dirs]

# create output directory
outdir = f"analysis/figures/{animal}/cell_noise_plots"
table_outdir = f"{outdir}/tables/"
isExist = os.path.exists(table_outdir)
if not isExist:
    print(f"creating output directory: {table_outdir}")
    os.makedirs(table_outdir)

# get the noise tables for each cell
for cell in cellnames:
    print(cell)
    files = [f"results/{animal}/{cell}/{cell}_0noiseF_Table.csv",f"results/{animal}/{cell}/{cell}_pinkF_Table.csv",f"results/{animal}/{cell}/{cell}_whiteF_Table.csv"]
    #create table file
    table_path = f"{table_outdir}/{cell}_noise.csv"
    with open(table_path,"w") as table:
        header_row = ""
        for file in files:
            if os.path.exists(file):
                noise = file.split("_")[1]
                if noise == "0noiseF":
                    noise_label="No Noise"
                if noise == "pinkF":
                    noise_label="Pink Noise"
                if noise == "whiteF":
                    noise_label="White Noise"
                with open(file,"r") as infile:
                    line = next(infile)
                    if not header_row:
                        header_row = f"noise,{line}"
                        table.write(header_row)
                    skip = next(infile) #empty row
                    line = next(infile)
                    line = f'{noise_label},{line}'
                    table.write(line)
    table.close()

    # read table and prepare data
    data = pd.read_csv(table_path)
    #if data frame is empty move to next cell
    if data.empty:
        continue 
    df = pd.DataFrame()

    treat=data.noise
    current=""
    for t in treat:
        subdata = data[data.noise == t]
        subdata.reset_index(drop=True,inplace=True)
        #get current injection
        if not current:
            CIstr=subdata["current_injection[][0]"].astype(str)
            CIstr1=CIstr[0].replace("[","").replace("]","")
            CIstr2=re.split(" ",CIstr1)
            CInum=[int(n) for n in CIstr2]
            df["CI"]=CInum
            current = CInum
        #get AP data
        APstr=subdata["APs[][0]"].astype(str)
        #remove brackets
        APstr1=APstr[0].replace("[","").replace("]","")
        #split numbers into list
        APstr2=re.split(" ",APstr1)
        #convert character into integer
        APnum=[int(n) for n in APstr2]
        #add values to dataframe
        df[t]=APnum
        
    #plot current injection vs APs
    df.plot(x="CI",marker="o")
    plt.title(cell)
    plt.xlabel("Current Injection")
    plt.ylabel("Number of APs")
    plt.xticks(list(range(100,1100,100)))
    plt.savefig(f'{outdir}/{cell}.png')
    plt.close()

os.system(f"touch {table_outdir}/done.txt")
