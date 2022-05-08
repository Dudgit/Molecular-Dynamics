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




def readXYZ(filename):
    itemsNa=[]
    itemsF=[]
    boxZ = []

    with open(filename, mode='r') as f:
        sor = f.readline()
        sor = f.readline()
        runnum = sor.split()[-1]

        runnum = int(runnum)

        #box data, nem fontos
        sor = f.readline()
        sor = f.readline()
        sor = f.readline()
        sor = f.readline()
        sor = f.readline()
        sor = f.readline()
        
        sor = sor.split()
        boxZ.append(float(sor[1])-float(sor[0]))

        cimkek = f.readline().split()
        cimkek = cimkek[2:]
        cimkek.append('timestep')

        dataF = f.readline().split()
        dataNa = f.readline().split()
        dataF.append(runnum)
        dataNa.append(runnum)

        itemsF.append(dataF)
        itemsNa.append(dataNa)

        while(f.readline()):
            runnum = f.readline().split()
            runnum = int(runnum[0])
            #print(runnum)
            sor = f.readline()
            sor = f.readline()
            sor = f.readline()
            sor = f.readline()
            sor = f.readline()
            sor = f.readline()
            sor = f.readline()
            dataF = f.readline().split()
            dataNa = f.readline().split()
            dataF.append(runnum)
            dataNa.append(runnum)
            itemsF.append(dataF)
            itemsNa.append(dataNa)
            
    frameNa = pd.DataFrame(data=itemsNa, columns=cimkek)
    frameF = pd.DataFrame(data=itemsF, columns=cimkek)
    frameNa = frameNa.apply(pd.to_numeric)
    frameF = frameF.apply(pd.to_numeric)
    
    frameF['lejjebbZ'] = frameF.z-boxZ
    frameF['feljebbZ'] = frameF.z+boxZ
    
    return(frameNa, frameF)
    
def proj(frame1, posDiff, sign='+'):
    prod = []
    for (x1, y1, z1, x2, y2, z2, norm) in zip(frame1.fx, frame1.fy, frame1.fz, 
                                              posDiff.x, posDiff.y, posDiff.z, posDiff['abs']):
        length = (x1*x2+y1*y2+z1*z2)/norm
        if sign=='-':
            prod.append(-length)
        else:
            prod.append(length)
    return prod


def shallOffset(posDiff, posDiff_le, posDiff_fel, ionDist):
    ans = []
    for (orig, le, _) in zip(posDiff, posDiff_le, posDiff_fel):
        if abs(orig-ionDist) < 0.5:
            ans.append('orig')
        elif abs(le-ionDist) <0.5:
            ans.append('le')
        else:
            ans.append('fel')
    return ans

def applyOffset(col_orig, col_le, col_fel, col_offset):
    if col_offset == 'le':
        true_val = col_le
    elif col_offset =='fel':
        true_val = col_fel        
    else:
        true_val = col_orig        
    return true_val

def ionForce(filename, ionDist):
    (Na,F) = readXYZ(filename)
    
    Na['f'] = (Na.fx**2 + Na.fy**2 + Na.fz**2)**0.5
    F['f'] = (F.fx**2 + F.fy**2 + F.fz**2)**0.5
    
    posDiff = pd.DataFrame(data=Na.x-F.x, columns=['x'])
    posDiff['y'] = Na.y-F.y
    posDiff['z_orig'] = Na.z-F.z
    posDiff['z_fel'] = Na.z-F.feljebbZ
    posDiff['z_le'] = Na.z-F.lejjebbZ
    posDiff['abs_orig'] = (posDiff.x**2+ posDiff.y**2 + posDiff.z_orig**2)**0.5
    posDiff['abs_fel'] = (posDiff.x**2+ posDiff.y**2 + posDiff.z_fel**2)**0.5
    posDiff['abs_le'] = (posDiff.x**2+ posDiff.y**2 + posDiff.z_le**2)**0.5
    
    posDiff['offset'] = shallOffset(posDiff.abs_orig, posDiff.abs_le, posDiff.abs_fel, ionDist)
    
    posDiff['abs'] = posDiff.apply(lambda x: applyOffset(x.abs_orig, x.abs_le, x.abs_fel, x.offset), axis=1)
    posDiff['z'] = posDiff.apply(lambda x: applyOffset(x.z_orig, x.z_le, x.z_fel, x.offset), axis=1)

    posDiff = posDiff.drop(['abs_orig', 'abs_fel', 'abs_le', 'z_orig', 'z_le', 'z_fel'], axis=1)
    
    Na['fproj_f'] = proj(Na, posDiff, sign='-')
    F['fproj_na'] = proj(F, posDiff)
    
    return(F['fproj_na'].mean(), Na['fproj_f'].mean())

