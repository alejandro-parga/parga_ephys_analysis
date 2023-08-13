import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import seaborn as sns

animal = snakemake.params.animal
prop = snakemake.params.prop
outfile = snakemake.output.out
outpath = os.path.dirname(outfile)

# hard coded variables for testing
# animal = "NHP"
# prop = "FITable"
# outfile = f"analysis/figures/{animal}/drug/{prop}/done.txt"
# outpath = os.path.dirname(outfile)

print(f"animal: {animal}")
print(f"table: {prop}")
print(f"outfile = {outfile}")
print(f"outpath = {outpath}")

# get list of file paths base on variables
files = [f"analysis/{animal}/XE991a_{prop}_merged.csv",f"analysis/{animal}/XE991b_{prop}_merged.csv"]

#check if output directory exists
#if not, create
isExist = os.path.exists(outpath)
if not isExist:
    print(f"creating output directory: {outpath}")
    os.makedirs(outpath)

#initialize empty data frame
df = pd.DataFrame()
for file in files:
    #read data, fills empty with NaN
    data = pd.read_csv(file)
    #add column for drug
    drug = re.split("XE991|_", file)[1]
    if(drug == "a"):
        drug_label = "baseline"
    if(drug == "b"):
        drug_label = "drug"
    data["drug"]=drug_label

    #for the FI table average/count for the columns with multiple values 
    if prop == "FITable":
        #add column for AP counts/mean
        apcount = list()
        apmean = list()
        for str in data["APs[][0]"]:
            str1 = str.replace("[","").replace("]","")
            str2 = re.split(" ",str1)
            num = [int(n) for n in str2]
            count = sum(num)
            avg = count/len(num)
            apcount.append(count)
            apmean.append(avg)

        #add column for AP Hz
        apHzmean = list()
        for str in data["maxHz[][0]"]:
            str1 = str.replace("[","").replace("]","")
            str2 = re.split(" ",str1)
            num = [float(n) for n in str2]
            count = sum(num)
            avg = count/len(num)
            apHzmean.append(avg)

        data["APcount"] = apcount
        data["APmean"] = apmean
        data["APmeanHz"] = apHzmean

    df = pd.concat([df,data])
    df.reset_index(drop=True,inplace=True)

#plot each property by drug
colnames = list(df.columns)

# plot each column by drug: 
skip = ['cell_name','APs[][0]', 'maxHz[][0]','current_injection[][0]','drug']
for c in colnames[1:]:
    if c in skip:
        continue
    
    print(f"plotting {c}")

    #check for missing values
    #add to verify list and skip
    if df[c].isnull().values.any():
        nulldf = df[df[c].isnull()]
        sub_nulldf = nulldf[["cell_name"]]
        sub_nulldf["table"] = prop
        sub_nulldf["metric"] = c

        sub_nulldf.to_csv("verify/check_drug.txt",index=False,mode="a")

    #plot!
    sns.set(style="whitegrid")
    p = sns.stripplot(x="drug", y=c, data=df,size=8)
    sns.boxplot(showfliers = False,                     
                showmeans=True,
                meanline=True,
                meanprops={'color': 'k', 'ls': '-', 'lw': 2},
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                #zorder=10,
                x="drug",
                y=c,
                data=df,
                showbox=False,
                showcaps=False,
                ax=p,
                width=0.5)
    plt.xlabel("")
    plt.title(c)
    plt.savefig(f'{outpath}/{c}.png')
    plt.close()

f=open(outfile, "w")
f.write("Done!")
