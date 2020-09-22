# nanoparticle_assembly

This is the project repository for the work titled : #Examining the Self-Assembly of Patchy Alkane-Grafted Silica Nanoparticles using Molecular Simulation

Here contains the procedures used to produce the data and figures in this paper. Cloning this notebook, along with the relevant raw data available at Zenodo:DOI, will allow one to replicate all of the data extracted from the hoomd-blue trajectories. Unfortunately, due to updates to the hoomd-blue codebase, the actual simulation workflow cannot be replicated with the code provided here. It is still provided in order to provide relvant inputs that went into running those simulations.

*INSTALLING AND RUNNING THESE NOTEBOOKs*
After cloning the repo, enter into the Setup Folder. Open the READAME.md file and follow the instructions to set up your Conda environment.

Important notebooks to look at:
- Plot_Data.ipynb is the notebook that plots all of the finalized values from the .csv tabulated metrics
- Read_Trajectory.ipynb reads in data from the trajectories, and creates center's of mass for running calculations.
- Tethered_Particles_Calculations.ipynb takes the data read in from the full trajectories, and does calculations for the phases
- vis_scripts.ipynb makes hoomdxml files that can be opened and visualized in VMD for seeing the full structures
- R_g_Asphericity_Calcs_for_full_system.ipynb calculate the radius of gyration and asphericity from the single nanoparticles.

