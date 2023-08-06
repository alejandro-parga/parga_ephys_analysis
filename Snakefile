# Snakefile
# created by Alejandro Parga
# updated July 2023
# compiling csv tables from MIES

### VARIABLES
#memprop_tables = ["FITable","VITable","ZAPTable","APTable"]
noise_tables = ["0noiseF_Table","0noiseF1_Table","0noiseF2_Table","pinkF_Table","pinkF1_Table","pinkF2_Table","whiteF_Table","whiteF1_Table","whiteF2_Table"]
animals = ["NHP","Thy1","Etv1"]
fs = ["F","F1","F2"]
drugs = ["a","b"]

rule all:
    input:
        expand("analysis/{animal}/{table}_merged.csv",animal = animals,table = ["FITable","VITable","ZAPTable","APTable"]),
        expand("analysis/figures/{animal}/{table}/done.txt",animal = animals,table = ["VITable","ZAPTable"]),
        expand("analysis/{animal}/{table}_merged.csv",animal = animals,table = noise_tables),
        expand("analysis/figures/{animal}/frequency_{f}.png",animal = ["NHP","Thy1"],f=fs),
        #expand("analysis/figures/{animal}/cell_noise_plots/{cell}_noise{f}.png",animal = ["NHP","Thy1"],f=fs)
        #fix this later
        expand("analysis/figures/{animal}/cell_noise_plots/tables/done{f}.txt",animal = ["NHP","Thy1"],f=fs),
        expand("analysis/{animal}/XE991{drug}_{table}_merged.csv", animal = ["NHP"], drug = drugs, table = ["FITable","VITable","ZAPTable"])

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
rule plot_noise_F_by_cell:
    input:
        # "results/{animal}/{cell}/{cell}_0noise{f}_Table.csv",
        # "results/{animal}/{cell}/{cell}_pink{f}_Table.csv",
        # "results/{animal}/{cell}/{cell}_white{f}_Table.csv"
        "results/{animal}/"
    output:
        #fix this later
        "analysis/figures/{animal}/cell_noise_plots/tables/done{f}.txt"
        #"analysis/figures/{animal}/cell_noise_plots/{cell}_noise{f}.png"
    params:
        animal="{animal}",
        f="{f}"
    wildcard_constraints:
        animal="NHP|Thy1"
    script:
        "scripts/plot_noise_per_cell.py"

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
# rule plot_drug:
#     input:
#         "analysis/{animal}/0noise{f}_Table_merged.csv",
#         "analysis/{animal}/pink{f}_Table_merged.csv",
#         "analysis/{animal}/white{f}_Table_merged.csv"
#     output:
#         "analysis/figures/{animal}/frequency_{f}.png"
#     params:
#         animal="{animal}",
#         f="{f}"
#     wildcard_constraints:
#         animal="NHP|Thy1"
#     script:
#         "scripts/plot_noise_tables.py"

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
