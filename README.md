# trophic-tuna
Molecular Dynamics simulation of the H88Y substitution mutation in cardiac actin (ACTC1), using GROMACS on a HPC cluster via the SLURM workload manager.
***
## Subheading ##
Nomenclature: "6DJN-1-H88Y" would mean the 6DJN pdb file with the H88Y mutation added, as a single subunit. "6DJO-8-WT" would mean the wildtype 6DJO pdb file as an 8-subunit-long filament. 
Note: The filament takes substantially more time to perform MD on, and in retrospect 8 subunits was overkill for a primary analysis. I ended up using a 6-subunit-long filament instead, but the commands are the same for the protomer as for the polymers. The only difference to keep in mind is that the parameters (and hence parameter files) differ, though I have included both files in the repository. 

The entire process from pre-processing to the MD simulation itself is split into 3 major steps:
1. **Pre-processing and energy minimization**. Here, the protein is centered in a cubic system where water molecules and ions are added to balance any net charge the protein may have. The protein is set up such that if any atoms hit one face of the cube, rather than bounce off (that would simulate protein-surface physics, which are interesting but don't represent most interactions in a cell) they pass through the face and exit the other side. The sides of the cube are far enough from the center that the protein will never interact with itself. The system energy is calculated and steadily decreased until it reaches a local minimum. One thing to be aware of is that for your system, that local minimum may not be the global minimum, and making the step-size too small increases that risk.

2. **Temperature and Pressure Equilibration**. After the starting structure of the protein is more representative of what we're likely to see in a cell rather than the easiest structure to crystallize (i.e. the pdb structure), the system must be brought from no kinetic/thermal energy up to atmospheric pressure and standard temperature (1 atm and 298 K). I was shown by a collegue that using a temperature of 310 Kelvin decreased simulation time without significant sacrifices to accuracy, which is why my parameter files set the reference temperature to 310 Kelvin. For initial velocities, I used pseudo-random sampling from a Maxwell-Boltzmann distribution.

3.  **Molecular Dynamics Simulation**. Finally the system is ready for MD simulation. I used time-steps of 2 fs and a total simulation time of 100 ns.

***

## Analyses ##

I performed a few basic analysis like RMSD, ΔRMSD/t, RMSF, radius of gyration, and principle components analysis. I also want to perform dynamic network analysis (for the allosteric network) but haven't made the time yet. Code for getting the above measurements is included in an analysis.sh file, but I did not include my R code because it was just graphing and many tutorials for that can be found online. I used the Peptides library (mainly readXVG()) and the FisherZ function from DescTools (Fisher-Transformation for Correlation to z-Score), both of which I found very handy.

***

I am happy to discuss any other details a user might want to know about the project, just send me an email :) noahzeidenbergATgmail.com
