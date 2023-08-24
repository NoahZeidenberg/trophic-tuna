#!/bin/bash
#SBATCH --time=0-00:45:00           # time limit (D-HH:MM:SS)
#SBATCH --account=[[YOUR ACCOUNT]]
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=0 
#SBATCH --output=6DJN_1_H88Y_pre_min.out   #6DJN (pdb), single (1) subunit, ACTC1-H88Y
#SBATCH --job-name=6DJN_1_H88Y_pre_min
#SBATCH --mail-user=[[YOUR EMAIL]]
#SBATCH --mail-type=FAIL            # receive email notifications if the script prompts an error 
#SBATCH --mail-type=ALL             # receive email notifications if the script runs successfully

module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 gromacs/2022.3
export GMXLIB=/home/[your_account]/gromacs-2022.3/share/top
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"   #optimized distribution of tasks based on available CPUs

#queue NPT + NVT equilibration(s) to run if energy minimization (this script) finishes successfully

sbatch -d afterok:$SLURM_JOB_ID --chdir=/home/[your_account]/scratch/equilibration 6DJN-1-H88Y_equilib.sh

#generate topology + position restraint

#the pdb file is in a folder called "pdb_files", the output files for 
# preprocessing/energy minimization are in a folder called "pre_min", 
# the topology files are in a folder called "topology_files", and the
# parameter files are in a folder called "parameter_files".

gmx pdb2gmx -f /home/[your_account]/scratch/pdb_files/6DJN-1-H88Y.pdb -o /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_processed.gro -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top -i /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_posre.itp -ff charmm36-jul2022 -water tip3p -merge all

#using the CHARMM36 forcefield as it is optimal for single-protein systems (i.e. rather than protein-DNA, protein-RNA, protein-lipid, etc. systems).
#using the TIP3P explicit water model as it is more accurate than an implicit model, and the HPC cluster can handle the higher computational load.

#define box around protein

gmx editconf -f /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_processed.gro -o /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_newbox.gro -c -d 1.0 -bt cubic

#defining a cubic box with faces at least 1 nm away from all protein atoms, with the protein at the center

#solvate system 

gmx solvate -cp /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_newbox.gro -cs spc216.gro -o /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_solv.gro -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top

#Add ions to system for net 0 charge on protein

gmx grompp -f /home/[your_account]/scratch/parameter_files/ions.mdp -c /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_solv.gro -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top -o /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_ions.tpr -maxwarn 3

echo SOL | gmx genion -s /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_ions.tpr -o /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_solv_ions.gro -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top -pname NA -nname CL -neutral

#                                                                            sodium for + charge, chlorine for - charge.

#energy minimization by steepest descent

gmx grompp -f /home/[your_account]/scratch/parameter_files/minim_sd.mdp -c /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_solv_ions.gro -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top -o /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_em.tpr -maxwarn 3

gmx mdrun -s /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_em.tpr -deffnm 6DJN-1-H88Y -c /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_em.gro

#create index file

echo -e '"Protein"|"MG"|"ADP"|"PO4" \n "Water"|"NA" \n q' | gmx make_ndx -f 6DJN-1-H88Y_em.gro -o 6DJN-1-H88Y_index.ndx