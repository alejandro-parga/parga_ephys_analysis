# master merged
# assign varibale animal
anmial = snakemake.params.animal


# read in FITable_merged.csv, VITable_merged.csv, ZAPTable_merged.csv, APTable_merged.csv
analysis/{animal}/FITable_merged.csv
# append the merged tables by cellname
# fill in N/A
# wrtie the output