import csv
import os
import glob

# input variables from snakefile
dir_path = snakemake.input.rec
prop = snakemake.params.prop
outpath = snakemake.output.merged

# #hardcoded variables for testing
# print("TESTING")
# prop = "FITable"
# dir_path = "results/Thy1/"
# outpath = f"analysis/Thy1/{prop}_merged.csv"

#for log, print variables
print(f"property: {prop}")
print(f"input data: {dir_path}")
print(f"output file: {outpath}")

#check if output directory exists
#if not, create
outdir = os.path.dirname(outpath)
isExist = os.path.exists(outdir)
if not isExist:
    print(f"creating output directory: {outdir}")
    os.makedirs(outdir)

# get the files for this animal/property
file_path = f"{dir_path}/*/*_{prop}.csv"
files = glob.glob(file_path)

# read files and 
# assumes all files have the same columns
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