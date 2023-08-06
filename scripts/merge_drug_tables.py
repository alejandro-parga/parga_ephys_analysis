import csv
import os
import glob

# hardcoded variables for testing
# print("TESTING")
# prop = "FITable"
# dir_path = "results/NHP"
# drug = "a"
# outpath = "analysis/NHP/XE991a_FITable_merged.csv"

# input variables from snakefile
dir_path = snakemake.input.rec
prop = snakemake.params.prop
drug = snakemake.params.drug
outpath = snakemake.output.merged

# for log, print variables
print(f"property: {prop}")
print(f"input data: {dir_path}")
print(f"drug: {drug}")
print(f"outpath file: {outpath}")

# check if output directory exists
# if not, create
outdir = os.path.dirname(outpath)
if not os.path.exists(outdir): 
    print(f"creating output directory: {outdir}")
    os.makedirs(outdir)

# getting list of files
filepath = f"{dir_path}/*/*_XE991{drug}_{prop}.csv"
files = glob.glob(filepath)

# read the files
with open(outpath, "w") as outfile:
    header_row = ""
    for f in files:
        print(f"reading file: {f}")
        with open(f, "r") as infile:
            line = next(infile) #throw out header row except first time
            if not header_row:
                header_row = line
                outfile.write(header_row)
            if prop == "FITable":
                line = next(infile)
            for line in infile:
                outfile.write(line)

print("done!")

