#!/bin/bash
#SBATCH --time=0-00:20:00           # time limit (D-HH:MM:SS)
#SBATCH --account=[[YOUR ACCOUNT]]
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=0 
#SBATCH --output=6DJN_1_H88Y_md_analysis.out
#SBATCH --job-name=6DJN-1-H88Y_md_analysis
#SBATCH --mail-user=[[YOUR EMAIL]]
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 gromacs/2022.3
export GMXLIB=/home/[your_account]/gromacs-2022.3/share/top
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

#process trajectory files:

#select protein atom subset and make molecule whole
echo -e 'Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_whole.xtc -pbc whole

#cluster atoms in protein index
echo -e 'Protein \n Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_whole.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_cluster.xtc -pbc cluster

#extract starting frame with initial velocities (.gro, not .pdb)
echo -e 'Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_cluster.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_cluster.gro -dump 0

#remove jumps across periodic boundary conditions (PBCs)
echo -e 'Protein \n Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_cluster.gro -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_cluster.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_nojump.xtc -center -pbc nojump

#center protein in the box
echo -e 'Protein \n Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_nojump.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_center.xtc -center -ur compact -pbc mol

#dump centered structure (.gro)
echo -e 'Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_center.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_center.gro -dump 0

#fit protein to dumped reference structure, update, and fit iteratively...
echo -e 'Protein \n Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_center.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.xtc -fit progressive

#dump first frame for visualization in vmd 
echo -e 'Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.xtc -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.gro -dump 0

#dump last frame for pymol (.pdb)
echo -e 'Protein \n' | gmx trjconv -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.gro -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.pdb -dump 50000000

#------------------------------------------------------------------------------------

#Get RMSD
echo -e 'Backbone \n' | gmx rms -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.xtc -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_rmsd.xvg -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx

#Get RMSF (averaged for each residue)
echo -e 'Protein \n' | gmx rmsf -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.xtc -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_rmsf.xvg -ox /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_rmsf.pdb -res -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx

#Radius of Gyration
echo -e 'Protein \n' | gmx gyrate -s /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_md_trajec_fit.xtc -o /home/[your_account]/scratch/MD_analysis/6DJN-1-H88Y_gyrad.xvg -p -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx

#Get first ten principle components of protein trajectory, can be modeled easily in PyMOL
echo -e 'Backbone \n Backbone \n' | gmx covar -s /home/[your_account]/scratch/MD/6DJO-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJO-1-H88Y_md_trajec_fit.xtc -o /home/[your_account]/scratch/MD_analysis/6DJO-1-H88Y_eigen.xvg -last 10

echo -e 'Backbone \n Backbone \n' | gmx anaeig -s /home/[your_account]/scratch/MD/6DJO-1-H88Y_md.tpr -f /home/[your_account]/scratch/MD_analysis/6DJO-1-H88Y_md_trajec_fit.xtc -extr /home/[your_account]/scratch/MD_analysis/6DJO-1-H88Y_extreme10.gro -first 1 -last 10 -nframes 50


