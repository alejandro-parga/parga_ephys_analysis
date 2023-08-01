import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

animal = snakemake.params.animal
f = snakemake.params.f

#hard coded variables for testing
# print("TESTING")
# animal = "Thy1" #options: NHP,Thy1
# f = "F1" #options: F,F1,F2

# get file list based on 
files = [f'analysis/{animal}/0noise{f}_Table_merged.csv',f'analysis/{animal}/pink{f}_Table_merged.csv',f'analysis/{animal}/white{f}_Table_merged.csv']

outpath = f"analysis/figures/{animal}/"

#check if output directory exists
#if not, create
isExist = os.path.exists(outpath)
if not isExist:
    print(f"creating output directory: {outpath}")
    os.makedirs(outpath)

#initialize data frame with x values
df = pd.DataFrame(list(range(100,1100,100)),columns=['x'])
for file in files:
    #read data, fills empty with NaN
    data = pd.read_csv(file)

    #get noise type from file name
    #noise = re.split("/|_",file)[2]
    #noise = re.sub("F","",noise) # remove F,F1,F2
    noise = file.split("/")[2]
    noise = noise.split("_")[0]

    #get just the y values
    y_data = data.filter(regex="y_")
    # calculate the mean for each y value, by default the function ingores NaN values
    y_mean = y_data.mean()
    # calculate standard error of the mean 
    # set ddof=1 to use sample standard deviation instead of population standard deviation
    # np.std(data,ddof=1) / np.sqrt(np.size(data))
    y_se = np.std(y_data,ddof=1)/np.sqrt(np.size(y_data))

    # add mean and standard error to dataframe
    df[f'{noise}_mean'] = list(y_mean)
    df[f'{noise}_se'] = list(y_se)


#plot means and standard errors
df.plot(x="x",y=[f"0noise{f}_mean",f"pink{f}_mean",f"white{f}_mean"],kind="line",marker='.')
plt.legend(loc=2)
plt.errorbar("x",f"0noise{f}_mean",yerr = f"0noise{f}_se",data=df,fmt=" ",color="k")
plt.errorbar("x",f"pink{f}_mean",yerr = f"pink{f}_se",data=df,fmt=" ",color="k")
plt.errorbar("x",f"white{f}_mean",yerr = f"white{f}_se",data=df,fmt=" ",color="k")
plt.title(animal)
plt.xlabel("")
plt.ylabel("Frequency (Hz)")
plt.xticks(list(range(100,1100,100)))
plt.savefig(f"{outpath}/frequency_{f}.png")
plt.close()