from __future__ import division

import mdtraj as md
import numpy as np

def accessible_surface_area(traj, nanoparticle_radius, probe_radius=0.25):
    radii_exchange = {'C': nanoparticle_radius, 'N': MMM_sigma/2, 'O': MME_sigma/2}

    sasa = md.shrake_rupley(traj, mode='residue', change_radii=radii_exchange,
                            probe_radius=probe_radius, n_sphere_points=1000)

    return np.column_stack((traj.time, sasa[:,0]))

