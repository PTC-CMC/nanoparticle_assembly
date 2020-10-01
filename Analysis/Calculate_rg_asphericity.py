import numpy as np
import pandas
import os
import scipy
import MDAnalysis as md

path = '../trajectories/single_particle'
dcdpath = 'tnp.dcd'
toppath = 'final.mol2'
statepath= 'signac_statepoint.json'


chain_rg=[]
full_rg=[]
asphere=[]

for filename in os.listdir(os.path.join(os.getcwd(),path)):
    full_dcdpath = os.path.join(os.path.join(path,filename),dcdpath)
    full_toppath = os.path.join(os.path.join(path,filename),toppath)
    full_statepath = os.path.join(os.path.join(path,filename),statepath)

    print(full_toppath + '\n\n\n')
    universe = md.Universe(full_toppath,full_dcdpath,topology_format='MOL2')
    atom_length = len(universe.atoms)
    chain_ids = "bynum 154:" + str(atom_length)
    full_ids="bynum 154:" + str(atom_length)
    with open(full_statepath,'r') as myrfile:
        lines = myrfile.readlines()
        chainlength = int(lines[4][17:19])
    atomg = universe.atoms.select_atoms(chain_ids)
    total = 0
    #for chain_id in np.arange(0,len(atomg)/chainlength*3,1):
    for chain_id in np.arange(0,len(atomg)/chainlength*3,1):
        chain = atomg[int(chain_id*(chainlength/3)):int((chain_id+1)*chainlength/3)]
        chain.positions=chain.positions*0.395
        total+=chain.radius_of_gyration(pbc=True)
        #print(len(chain))    
    chain_rg.append(total/(len(atomg)/chainlength*3))                   
    full_rg.append(universe.atoms.radius_of_gyration()*0.395)
    asphere.append(universe.atoms.asphericity())
df = pandas.DataFrame({'Graft_rg': 2.5/np.array(chain_rg),
                       'NP_rg': 2.5/np.array(full_rg),
                       'Asphericity': asphere})
print(universe.atoms.select_atoms('bynum 155').total_mass(), len(atomg),chainlength,atomg.total_mass()/len(atomg))
df.to_csv('rg_asphere.csv')
