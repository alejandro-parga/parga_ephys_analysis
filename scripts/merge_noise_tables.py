import os
import glob

# input variables from snakefile
dir_path = snakemake.input.rec
prop = snakemake.params.prop
outpath = snakemake.output.merged

# #hardcoded variables for testing
# print("TESTING")
# prop = "0noiseF_Table"
# dir_path = "results/NHP/"
# outpath = f"analysis/NHP/{prop}_merged.csv"

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

# read files, skip empty line, and merge
with open(outpath, "w") as outfile:
    header_row = ""
    for f in files:
        print(f"reading file: {f}")
        # cell_name = f.split("/")[2] #in case I need F1, F2
        with open(f,"r") as infile:
            line = next(infile) #throw out header row except first time
            if not header_row:
                #create header row for the max number of columns and add cell name
                #if I need F1, F2 line below
                # header_row = 'cell_name,x_100,y_100,x_200,y_200,x_300,y_300,x_400,y_400,x_500,y_500,x_600,y_600,x_700,y_700,x_800,y_800,x_900,y_900,x_1000,y_1000\n'
                header_row = line
                outfile.write(header_row)
            empty_line = next(infile) #skip empty line
            line = next(infile)
            # line = f'{cell_name},{line}' #add cell name variable for F1 and F2
            outfile.write(line)

                    