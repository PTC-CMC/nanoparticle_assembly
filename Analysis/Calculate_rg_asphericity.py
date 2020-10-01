import numpy as np
import pandas
import os
import scipy
import MDAnalysis as md

path = '../trajectories/single_particles/'
dcdpath = 'tnp.dcd'
toppath = 'final.mol2'


r_g=[]
asphere=[]

for filename in os.listdir(os.path.join(os.getcwd(),path)):
    full_dcdpath = os.path.join(os.path.join(path,filename),dcdpath)
    full_toppath = os.path.join(os.path.join(path,filename),toppath)

    print(full_toppath + '\n\n\n')
    universe = md.Universe(full_toppath,full_dcdpath,topology_format='MOL2')
    r_g.append(universe.atoms.radius_of_gyration())
    asphere.append(universe.atoms.asphericity())
#print(r_g,asphere)
df = pandas.DataFrame({'Radius of Gyration (nm)': r_g,
                       'Asphericity': asphere})

df.to_csv('rg_asphere.csv')
