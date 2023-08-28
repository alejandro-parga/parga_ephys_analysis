# Snakefile
# created by Alejandro Parga
# updated July 2023
# compiling csv tables from MIES

### VARIABLES
#memprop_tables = ["FITable","VITable","ZAPTable","APTable"]
# below variables for F1 and F2 merged noise tables
# noise_tables = ["0noiseF_Table","0noiseF1_Table","0noiseF2_Table","pinkF_Table","pinkF1_Table","pinkF2_Table","whiteF_Table","whiteF1_Table","whiteF2_Table"]
noise_tables = ['0noiseF_Table','pinkF_Table','whiteF_Table'] # noise variable for only F merged noise tables
animals = ["NHP","Thy1","Etv1"]
fs = ["F","F1","F2"]
drugs = ["a","b"]

rule all:
    input:
        #makes membraine property tables "merge_tables"
        expand("analysis/{animal}/{table}_merged.csv",animal = animals,table = ["FITable","VITable","ZAPTable","APTable"]),
        #plot membraine property "plot_membrane_properties"
        expand("analysis/figures/{animal}/{table}/done.txt",animal = animals,table = ["VITable","ZAPTable"]),
        #makes merged noise tables
        expand("analysis/{animal}/{table}_merged.csv",animal = ["NHP", "Thy1"],table = noise_tables),
        # makes plots for F1 and F2 noise
        # expand("analysis/figures/{animal}/frequency_{f}.png",animal = ["NHP","Thy1"],f=fs),
        # expand("analysis/figures/{animal}/cell_noise_plots/{cell}_noise{f}.png",animal = ["NHP","Thy1"],f=fs)
        #fix this later
        # plot noise by cell 
        # expand("analysis/figures/{animal}/cell_noise_plots/tables/done.txt",animal = ["NHP","Thy1"]),
        expand("analysis/{animal}/XE991{drug}_{table}_merged.csv", animal = ["NHP"], drug = drugs, table = ["FITable","VITable","ZAPTable"]),
        expand("analysis/figures/{animal}/drug/{table}/done.txt",animal=["NHP"],table=["FITable","VITable","ZAPTable"]),
        # plot total FI by CI
        # expand("analysis/figures/{animal}/FITable/done.txt",animal=animals),
        # plot FI noise per cell with stdCI
        expand("analysis/figures/{animal}/stdCI_noise_plots/tables/done.txt", animal=["NHP", "Thy1"])

# works for memprop_tables
rule merge_tables:
    input:
        rec="results/{animal}/"
    output:
        merged="analysis/{animal}/{table}_merged.csv"
    params:
        prop="{table}"
    wildcard_constraints:
        table="VITable|ZAPTable|FITable|APTable"
    log:
        "logs/merge_tables/{animal}_{table}.log"
    script:
        "scripts/merge_tables.py"

# works for ["VITable","ZAPTable"]
rule plot_membrane_properties:
    input:
        merged=rules.merge_tables.output.merged
    output:
        out="analysis/figures/{animal}/{table}/done.txt"
    wildcard_constraints:
        table="VITable|ZAPTable"
    log:
        "logs/plot_membrane_properties/{animal}_{table}.log"
    script:
        "scripts/plot_properties.py"

# works for noise_tables
rule merge_noise_table:
    input:
        rec="results/{animal}/"
    output:
        merged="analysis/{animal}/{table}_merged.csv"
    params:
        prop="{table}"
    wildcard_constraints:
        table="0noiseF_Table|0noiseF1_Table|0noiseF2_Table|pinkF_Table|pinkF1_Table|pinkF2_Table|whiteF_Table|whiteF1_Table|whiteF2_Table"
    log:
        "logs/merge_tables/{animal}_{table}.log"
    script:
        "scripts/merge_noise_tables.py"

# works for NHP and Thy1
rule plot_noise_F:
    input:
        "analysis/{animal}/0noise{f}_Table_merged.csv",
        "analysis/{animal}/pink{f}_Table_merged.csv",
        "analysis/{animal}/white{f}_Table_merged.csv"
    output:
        "analysis/figures/{animal}/frequency_{f}.png"
    params:
        animal="{animal}",
        f="{f}"
    wildcard_constraints:
        animal="NHP|Thy1"
    script:
        "scripts/plot_noise_tables.py"

# works for NHP and Thy1 by cell
rule plot_noise_by_cell:
    input:
        # "results/{animal}/{cell}/{cell}_0noiseF_Table.csv",
        # "results/{animal}/{cell}/{cell}_pinkF_Table.csv",
        # "results/{animal}/{cell}/{cell}_whiteF_Table.csv"
        "results/{animal}/"
    output:
        #fix this later
        "analysis/figures/{animal}/cell_noise_plots/tables/done.txt"
        #"analysis/figures/{animal}/cell_noise_plots/{cell}_noise.png"
    params:
        animal="{animal}"
    wildcard_constraints:
        animal="NHP|Thy1"
    script:
        "scripts/plot_noise_per_cell.py"

# works for FI noise
# plots APs vs current injection per cell
# current injection standarize (100-1000 pA)
rule plot_noise_stdCI_per_cell:
    input:
        "results/{animal}/"
    output:
        "analysis/figures/{animal}/stdCI_noise_plots/tables/done.txt"
    params:
        animal="{animal}"
    script:
        "scripts/plot_noise_stdCI_per_cell.py"
        
# works for drug tables
rule merge_drug_tables:
    input:
        rec="results/{animal}"
    output:
        merged="analysis/{animal}/XE991{drug}_{table}_merged.csv"
    params:
        prop="{table}",
        drug="{drug}"  
    wildcard_constraints:
        table="FITable|VITable|ZAPTable|APTable"
    log:
        "logs/merge_drug_tables/{animal}_{drug}_{table}.log"
    script:
        "scripts/merge_drug_tables.py"

# works for NHP drug applications
rule plot_drug:
    input:
        a="analysis/{animal}/XE991a_{table}_merged.csv",
        b="analysis/{animal}/XE991b_{table}_merged.csv"
    output:
        out="analysis/figures/{animal}/drug/{table}/done.txt"
    params:
        animal="{animal}",
        prop="{table}"
    log: 
        "logs/plot_drug_tables/{animal}_{table}.log"
    script:
        "scripts/plot_drug_tables.py"

# works for ploting APs from FITable
rule plot_APs:
    input: 
        "analysis/{animal}/FITable_merged.csv"
    output:
        out="analysis/figures/{animal}/FITable/done.txt"
    params:
        animal="{animal}"
    log:
        "logs/plot_APs/{animal}.log"
    script:
        "scripts/plot_FI_APs_per_cell.py"

rule update_rulegraph:
    input:
        "Snakefile"
    output:
        "rulegraph.png"
    shell:
        "snakemake --forceall --rulegraph | dot -Tpng > rulegraph.png"


# delete tables to force re-creation for update 
rule clean:
    input:
        "analysis/NHP",
        "analysis/Thy1",
        "analysis/Etv1"
    shell:
        '''
        rm -rf {input}
        '''
