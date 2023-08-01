import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

# hard coded variables for testing
# print ("TESTINGS")
# animal = "Thy1" #options:NHP, Thy1
# f = "F" #options: F,F1,F2

animal = snakemake.params.animal
f = snakemake.params.f

# get list of cells
cellnames = os.listdir(f"results/{animal}/")

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
    files = [f"results/{animal}/{cell}/{cell}_0noise{f}_Table.csv",f"results/{animal}/{cell}/{cell}_pink{f}_Table.csv",f"results/{animal}/{cell}/{cell}_white{f}_Table.csv"]
    #create table file
    table_path = f"{table_outdir}/{cell}_noise.csv"
    with open(table_path,"w") as table:
        header_row = 'noise,x_100,100,x_200,200,x_300,300,x_400,400,x_500,500,x_600,600,x_700,700,x_800,800,x_900,900,x_1000,1000\n'
        table.write(header_row)
        for file in files:
            if os.path.exists(file):
                noise = file.split("_")[1]
                with open(file,"r") as infile:
                    skip = next(infile) #header row
                    skip = next(infile) #empty row
                    line = next(infile)
                    line = f'{noise},{line}'
                    table.write(line)
    table.close()

    # read table and prepare data
    data = pd.read_csv(table_path)
    #if data frame is empty move to next cell
    if data.empty:
        continue    
    y_data = data.filter(regex="noise|^[0-9]+")
    df = y_data.T
    df.rename(columns=df.iloc[0],inplace=True)
    df.drop(df.index[0],inplace=True)

    #plot data
    df.plot(marker=".")
    plt.title(cell)
    plt.xlabel("")
    plt.ylabel("Frequency (Hz)")
    plt.savefig(f"{outdir}/{cell}_noise{f}.png")
    plt.close()

os.system(f"touch {table_outdir}/done{f}.txt")
