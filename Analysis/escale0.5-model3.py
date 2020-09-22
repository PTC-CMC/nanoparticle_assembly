from __future__ import division

import itertools

from math import pi, sqrt
import random
import sys
from hoomd_script import *
import numpy as np

def _gen_velocities(system):
    for p in system.particles:
        mass = p.mass
        sigma = sqrt(reduced_temp / mass)
        p.velocity = (random.gauss(0, sigma), random.gauss(0, sigma),
            random.gauss(0, sigma))

'''
##### RUN INFO #####
--------------------
'''
name = 'escale0.5-model3'

# Timesteps (femtoseconds)
ts_minimize = 0.01 # femtoseconds
ts_eq1 = 1.0
ts_eq2 = 4.0 # femtoseconds

# Run times (femtoseconds)
time_eq1 = 1e5
time_eq2 = 10e6

# Temperature (Kelvin)
temp = 300

# Initialize random number generator
seed = 12345
random.seed(seed)

epsilon_factor = 0.5

'''
##### Cutoffs #####
-------------------
'''
# Cutoff for chain-chain interactions (nm)
cutoff_chains = 1.4

# Cutoff for nanoparticle core interactions (nm)
cutoff_np = 5.0

# Cutoff for chain-nanoparticle cross interactions (nm)
cutoff_cross = 5.0

'''
##### Units #####
-----------------
'''
# Reference distance unit (nm)
ref_dist = 0.395

# Reference energy unit (kcal/mol)
ref_energy = 0.091493

# Reference mass unit (amu)
ref_mass = 1.0

# Reference time unit - derived
ref_time = 638.42

# Reference inverse temperature unit - derived
ref_kb = 0.0217198

# Reference force unit - derived
ref_force = ref_energy / ref_dist

# Reduced temperature
reduced_temp = ref_kb * temp

'''
##### Load the initial topology #####
-------------------------------------
'''
system = init.read_xml('C18-C18.hoomdxml', wrap_coordinates=True)

'''
##### Group definitions #####
-----------------------------
'''
all = group.all()
rigid = group.rigid()
nonrigid = group.nonrigid()

'''
##### Initialize velocities #####
---------------------------------
'''
_gen_velocities(system)

'''
##### Define neighborlists #####
--------------------------------
'''
# Define neighborlist for chain-chain interactions
nl_chains = nlist.cell()
nl_chains.set_params(r_buff=0.25/ref_dist, check_period=2)
nl_chains.reset_exclusions(exclusions=['bond', 'angle', 'dihedral', 'body'])

# Define neighborlist for nanoparticle core interactions
nl_np = nlist.cell()
nl_np.set_params(r_buff=0.25/ref_dist, check_period=2)
nl_np.reset_exclusions(exclusions=['bond', 'angle', 'dihedral', 'body'])

# Define neighborlist for chain-nanoparticle cross interactions
nl_np_chain = nlist.cell()
nl_np_chain.set_params(r_buff=0.25/ref_dist, check_period=2)
nl_np_chain.reset_exclusions(exclusions=['bond', 'angle', 'dihedral', 'body'])

'''
##### Define forcefield parameters #####
----------------------------------------
'''
def potential_cutoff(n, m, sigma):
    return ((n/m)*sigma**(n-m))**(1/(n-m))

mie_chains = pair.mie(r_cut=cutoff_chains/ref_dist, nlist=nl_chains)
mie_np = pair.mie(r_cut=cutoff_np/ref_dist, nlist=nl_np)
mie_np_chain = pair.mie(r_cut=cutoff_cross/ref_dist, nlist=nl_np_chain)

chain_particle_types = ['_MMM', '_MME', '_MMMA', '_MMEA']
chain1_particle_types = ['_MMM', '_MME']
chain2_particle_types = ['_MMMA', '_MMEA']
particle_types = ['_CGN', '_MMM', '_MME', '_MMMA', '_MMEA']

epsilons = {'_MMM-_MMM': 0.389092, '_MME-_MME': 0.434996, '_MME-_MMM': 0.411404, '_CGN-_MMM': 0.6360,
            '_CGN-_MME': 0.6788, '_CGN-_CGN': 0.9286}
sigmas = {'_MMM-_MMM': 0.4582, '_MME-_MME': 0.4662, '_MME-_MMM': 0.4662, '_CGN-_MMM': 0.5291, 
        '_CGN-_MME': 0.5331, '_CGN-_CGN': 0.6}
ns = {'_CGN-_MMM': 29.6574, '_CGN-_MME': 28.0821, '_MMM-_MMM': 9, '_MME-_MME': 9, '_MME-_MMM': 9,
      '_CGN-_CGN': 20.0087}
ms = {'_CGN-_MMM': 5.5835, '_CGN-_MME': 5.6179, '_MMM-_MMM': 6, '_MME-_MME': 6, '_MME-_MMM': 6,
      '_CGN-_CGN': 4.7578}

epsilons['_MMMA-_MMMA'] = epsilons['_MMM-_MMM'] * epsilon_factor
epsilons['_MMEA-_MMEA'] = epsilons['_MME-_MME'] * epsilon_factor
epsilons['_MMEA-_MMMA'] = epsilons['_MME-_MMM'] * epsilon_factor
sigmas['_MMMA-_MMMA'] = sigmas['_MMM-_MMM']
sigmas['_MMEA-_MMEA'] = sigmas['_MME-_MME']
sigmas['_MMEA-_MMMA'] = sigmas['_MME-_MMM']
ns['_MMMA-_MMMA'] = ns['_MMM-_MMM']
ns['_MMEA-_MMEA'] = ns['_MME-_MME']
ns['_MMEA-_MMMA'] = ns['_MME-_MMM']
ms['_MMMA-_MMMA'] = ms['_MMM-_MMM']
ms['_MMEA-_MMEA'] = ms['_MME-_MME']
ms['_MMEA-_MMMA'] = ms['_MME-_MMM']
epsilons['_CGN-_MMMA'] = epsilons['_CGN-_MMM'] * epsilon_factor
epsilons['_CGN-_MMEA'] = epsilons['_CGN-_MME'] * epsilon_factor
sigmas['_CGN-_MMMA'] = sigmas['_CGN-_MMM']
sigmas['_CGN-_MMEA'] = sigmas['_CGN-_MME']
ns['_CGN-_MMMA'] = ns['_CGN-_MMM']
ns['_CGN-_MMEA'] = ns['_CGN-_MME']
ms['_CGN-_MMMA'] = ms['_CGN-_MMM']
ms['_CGN-_MMEA'] = ms['_CGN-_MME']

epsilons['_MMM-_MMMA'] = (epsilons['_MMM-_MMM'] * epsilons['_MMMA-_MMMA']) ** 0.5
epsilons['_MME-_MMEA'] = (epsilons['_MME-_MME'] * epsilons['_MMEA-_MMEA']) ** 0.5
epsilons['_MMEA-_MMM'] = (epsilons['_MME-_MMM'] * epsilons['_MMEA-_MMMA']) ** 0.5
epsilons['_MME-_MMMA'] = (epsilons['_MME-_MMM'] * epsilons['_MMEA-_MMMA']) ** 0.5
sigmas['_MMM-_MMMA'] = (sigmas['_MMM-_MMM'] + sigmas['_MMMA-_MMMA']) / 2
sigmas['_MME-_MMEA'] = (sigmas['_MME-_MME'] + sigmas['_MMEA-_MMEA']) / 2
sigmas['_MMEA-_MMM'] = (sigmas['_MME-_MMM'] + sigmas['_MMEA-_MMMA']) / 2
sigmas['_MME-_MMMA'] = (sigmas['_MME-_MMM'] + sigmas['_MMEA-_MMMA']) / 2
ns['_MMM-_MMMA'] = ns['_MMM-_MMM']
ns['_MME-_MMEA'] = ns['_MME-_MME']
ns['_MMEA-_MMM'] = ns['_MME-_MMM']
ns['_MME-_MMMA'] = ns['_MME-_MMM']
ms['_MMM-_MMMA'] = ms['_MMM-_MMM']
ms['_MME-_MMEA'] = ms['_MME-_MME']
ms['_MMEA-_MMM'] = ms['_MME-_MMM']
ms['_MME-_MMMA'] = ms['_MME-_MMM']

for type1, type2 in itertools.combinations_with_replacement(particle_types, r=2):

    type1, type2 = sorted((type1, type2))
    print(type1, type2)

    # Set chain-chain interactions
    if type1 in chain_particle_types and type2 in chain_particle_types:
        if (type1 in chain1_particle_types and type2 in chain1_particle_types):
            mie_chains.pair_coeff.set(type1, type2, epsilon=epsilons[type1 + '-' + type2] / ref_energy,
                                      sigma=sigmas[type1 + '-' + type2] / ref_dist, n=ns[type1 + '-' + type2],
                                      m=ms[type1 + '-' + type2],
                                      r_cut=potential_cutoff(ns[type1 + '-' + type2], ms[type1 + '-' + type2],
                                      sigmas[type1 + '-' + type2]/ref_dist)) 
        else:
            mie_chains.pair_coeff.set(type1, type2, epsilon=epsilons[type1 + '-' + type2] / ref_energy,
                                      sigma=sigmas[type1 + '-' + type2] / ref_dist, n=ns[type1 + '-' + type2],
                                      m=ms[type1 + '-' + type2])

    else:
        mie_chains.pair_coeff.set(type1, type2, epsilon=0.0, sigma=0.0, n=ns[type1 + '-' + type2],
                                  m=ms[type1 + '-' + type2], r_cut=False)


    if type1 not in chain_particle_types and type2 in chain_particle_types:
        if type2 in chain1_particle_types:
             mie_np_chain.pair_coeff.set(type1, type2, epsilon=epsilons[type1 + '-' + type2] / ref_energy,
                                         sigma=sigmas[type1 + '-' + type2] / ref_dist, n=ns[type1 + '-' + type2],
                                         m=ms[type1 + '-' + type2], r_cut=potential_cutoff(ns[type1 + '-' + type2],
                                         ms[type1 + '-' + type2], sigmas[type1 + '-' + type2]/ref_dist))
        else:
            mie_np_chain.pair_coeff.set(type1, type2, epsilon=epsilons[type1 + '-' + type2] / ref_energy,
                                        sigma=sigmas[type1 + '-' + type2] / ref_dist, n=ns[type1 + '-' + type2],
                                        m=ms[type1 + '-' + type2])
    else:
        mie_np_chain.pair_coeff.set(type1, type2, epsilon=0.0, sigma=0.0, n=ns[type1 + '-' + type2],
                                    m=ms[type1 + '-' + type2], r_cut=False)

    if type1 not in chain_particle_types and type2 not in chain_particle_types:
        mie_np.pair_coeff.set(type1, type2, epsilon=epsilons[type1 + '-' + type2] / ref_energy,
                              sigma=sigmas[type1 + '-' + type2] / ref_dist, n=ns[type1 + '-' + type2],
                              m=ms[type1 + '-' + type2])
    else:
        mie_np.pair_coeff.set(type1, type2, epsilon=0.0, sigma=0.0, n=ns[type1 + '-' + type2],
                              m=ms[type1 + '-' + type2], r_cut=False)

harmonic_bond = bond.harmonic()
harmonic_bond.bond_coeff.set('_MMM-_MMM', k=1232.1/(ref_energy/(ref_dist**2)), r0=0.364/ref_dist)
harmonic_bond.bond_coeff.set('_MME-_MMM', k=1232.1/(ref_energy/(ref_dist**2)), r0=0.365/ref_dist)
harmonic_bond.bond_coeff.set('_MMMA-_MMMA', k=1232.1/(ref_energy/(ref_dist**2)), r0=0.364/ref_dist)
harmonic_bond.bond_coeff.set('_MMEA-_MMMA', k=1232.1/(ref_energy/(ref_dist**2)), r0=0.365/ref_dist)

harmonic_angle = angle.harmonic()
harmonic_angle.set_coeff('_MMM-_MMM-_MMM', k=2.3846/ref_energy, t0=3.01942)
harmonic_angle.set_coeff('_MME-_MMM-_MMM', k=2.3846/ref_energy, t0=3.05433)
harmonic_angle.set_coeff('_MMMA-_MMMA-_MMMA', k=2.3846/ref_energy, t0=3.01942)
harmonic_angle.set_coeff('_MMEA-_MMMA-_MMMA', k=2.3846/ref_energy, t0=3.05433)
'''
##### Perform an energy minimization #####
------------------------------------------
'''
minimize = integrate.mode_minimize_fire(group=nonrigid, dt=ts_minimize/ref_time)
while not(minimize.has_converged()):
    run(1e4)
'''
##### Resize the box #####
----------------------------------------
'''
if round(system.box.Lx, 4) != 101.2658:
    resize_variant = variant.linear_interp(points=[(0, system.box.Lx),
                                                   (time_eq1/ts_eq1, 101.2658)])
    update.box_resize(Lx=resize_variant, Ly=resize_variant, Lz=resize_variant)
'''
##### Equilibrate the system #####
----------------------------------
'''
integrate.mode_standard(dt=ts_eq1 / ref_time)
# NOTE: Need to tune the drag coefficient
integrator = integrate.langevin(group=nonrigid, T=reduced_temp, seed=seed)
run(time_eq1 / ts_eq1)

#integrator = integrate.nve(group=nonrigid)
compute.thermo(group=all)
analyzer = analyze.log(quantities=['temperature', 'potential_energy',
    'kinetic_energy'], period=1000, filename='{}.thermo'.format(name), overwrite=True)
traj = dump.dcd(filename='{}.dcd'.format(name), period=10000, overwrite=True)
integrate.mode_standard(dt=ts_eq2 / ref_time)

for temperature in [1000, 800, 600, 400]:
    reduced_temp = ref_kb * temperature
    nvt_rigid = integrate.nvt_rigid(group=rigid, T=reduced_temp, tau=100/ref_time)
    run(time_eq2 / ts_eq2)
    nvt_rigid.disable()

reduced_temp = ref_kb * 300 
nvt_rigid = integrate.nvt_rigid(group=rigid, T=reduced_temp, tau=100/ref_time)
run(50e6 / ts_eq2)

