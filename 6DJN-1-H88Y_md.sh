#!/bin/bash
#SBATCH --time=0-10:58:00  # D-HH:MM:SS
#SBATCH --account=[[YOUR ACCOUNT]]
#SBATCH --nodes=1
#SBATCH --gres=gpu:p100:2
#SBATCH --ntasks-per-node=4  #number of threads
#SBATCH --cpus-per-task=8    #8 PMI processes per node
#SBATCH --mem=0  #mem per node
#SBATCH --job-name=6DJN_1-H88Y-md
#SBATCH --output=6DJN_1_H88Y_md.out
#SBATCH --mail-user=[[YOUR EMAIL]]
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=ALL


module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 gromacs/2021.4
export GMXLIB=/home/[your_account]/gromacs-2022.3/share/top
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

#The MD sim won't finish in 11 hours. Since SLURM waits to run scripts that request more time, it makes more sense to use the checkpoints allowed in mdrun and split up the simulation into 11-hour chunks. This is a short script that runs MD again from where it left off, if it hasn't reached target nsteps yet. It may throw an error if there is no checkpoint yet, but you can comment it out for the first run then next time you submit the script remove the #'s.

# First, check the status of the simulation
gmx_output=$(gmx check -f 6DJN-1-H88Y_md_checkpoint.cpt 2>&1)

# Parse the output to determine whether the simulation has encountered an error
if [[ "$gmx_output" == *"Fatal error:"* ]]; then
  # If the simulation has encountered a fatal error, do nothing
  :
else
  # If the simulation has not completed (checked with the parameter file), submit the job again
  sbatch -d afternotok:$SLURM_JOB_ID 6DJN-1-H88Y_md.sh
fi

gmx grompp -f /home/[your_account]/scratch/parameter_files/6DJN-md.mdp -c /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_npt.gro -t /home/[your_account]/scratch/equilibration/6DJN-1-H88Y_npt.cpt -p /home/[your_account]/scratch/topology_files/6DJN-1-H88Y_topol.top -n /home/[your_account]/scratch/W23/6DJN/1/H88Y/pre_min/6DJN-1-H88Y_index.ndx -o /home/[your_account]/scratch/MD/6DJN-1-H88Y_md.tpr

#run MD until completion or timeout

mpiexec gmx_mpi mdrun -deffnm 6DJN-1-H88Y_md -cpi /home/[your_account]/scratch/MD/6DJN-1-H88Y_md_checkpoint.cpt -cpo /home/[your_account]/scratch/MD/6DJN-1-H88Y_md_checkpoint.cpt
