# This is the project repository for the work titled : Examining the Self-Assembly of Patchy Alkane-Grafted Silica Nanoparticles using Molecular Simulation

Here contains the procedures used to produce the data and figures in this paper. Cloning this notebook, along with the relevant raw data available at Zenodo: 10.5281/zenodo.4042749, will allow one to replicate all of the data extracted from the HOOMD-Blue trajectories. Unfortunately, due to updates to the HOOMD-Blue codebase, the actual simulation workflow cannot be replicated with the code provided here. It is still provided in order to provide relevant inputs that went into running those simulations. The simulations were originally run in HOOMD (1.3.1).

*INSTALLING AND RUNNING THESE NOTEBOOKs*
Go to the section titled **Downloading and Working with Datasets** in the SI for this paper to follow the instructions to set up your Conda environment.

Important notebooks to look at:
- Plot_Data.ipynb is the notebook that plots all of the finalized values from the .csv tabulated metrics
- Read_Trajectory.ipynb reads in data from the trajectories, and creates center's of mass for running calculations.
- COG_Single_NP's.ipynb creates the center of geometries of the full nanoparticles to be utilized for nearest neighbor and RDF calculations.
- vis_scripts.ipynb makes hoomdxml files that can be opened and visualized in VMD for seeing the full structures
- Calculate_rg_asphericity.py calculate the radius of gyration and asphericity from the single nanoparticles.

*note that no code based on trajectory information will run if the bulk and single_particle data from our Zenodo DOI is not downloaded and moved into the trajectory folder. This includes the Read_Trajectory.ipynb, COG_Single_NP's.ipynb, Read_Trajectory.ipynb, and vis_scripts.ipynb. Data analysis through Plot_Data.ipynb will be the only thing that runs fully just based off of data in this repository*
