#!/bin/bash
#SBATCH --time=1-23:58:00           # time limit (D-HH:MM:SS)
#SBATCH --account=[[YOUR ACCOUNT]]
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=0 
#SBATCH --output=6DJN_1_H88Y_equilib.out
#SBATCH --job-name=6DJN_1_H88Y_equilib
#SBATCH --mail-user=[[YOUR EMAIL]]
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 gromacs/2022.3
export GMXLIB=/home/[your_account]/gromacs-2022.3/share/top
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

#queue MD to run if equilibration finishes succesfully

sbatch -d afterok:$SLURM_JOB_ID --chdir=/home/[your_account]/scratch/MD 6DJN-1-H88Y_md.sh

#run temperature equilibrium

gmx grompp -f /home/[your_account]/scratch/parameter_files/6DJN-nvt.mdp -c /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_em.gro -r /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_em.gro -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_nvt.tpr

#use mpiexec and gmx_mpi for parallel computing
mpiexec gmx_mpi mdrun -deffnm /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_nvt

#run pressure equilibrium

gmx grompp -f /home/[your_account]/scratch/parameter_files/6DJN-npt.mdp -c /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_nvt.gro -r /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_nvt.gro -t /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_nvt.cpt -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top -n /home/[your_account]/scratch/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_npt.tpr

mpiexec gmx_mpi mdrun -deffnm /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_npt 
