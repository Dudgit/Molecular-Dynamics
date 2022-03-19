#%%
dataPath = "../data/"
from lammps import lammps
import pandas as pd
import matplotlib.pyplot as plt
import  plotsettings

def readRDF(fpath:str):
    """
    TODO: Description needed
    """
    return pd.read_csv(f"{dataPath}rdf.dat",skiprows=4, names=['binNum','binEdge','rdf','coordNum'],sep=" ",index_col= 0)



#plt.plot(df['binEdge'], df['rdf'], '-')
#plt.xlabel('bin edge [A]', fontsize=15)
#plt.ylabel('rdf', fontsize=15)
#plt.title('Radial Distribution Function of Oxygens (special_bond excluded)', fontsize=18)
