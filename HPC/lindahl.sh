#!/bin/bash
#SBATCH --time=1-23:59:59
#SBATCH --nodes=1
#SBATCH -p lindahl3,lindahl4
#SBATCH --job-name=yourName
#SBATCH --mail-user=yourEmail
#SBATCH --mail-type=ALL
#SBATCH -G 1

# LOAD MODULES

module load cmake/latest
module load gcc/7.4
module load cuda/10.2

simdir="$PWD"

# COMPILE CUSTOM GROMACS ON SCRATCH

builddir="/scratch/$USER/$SLURM_JOBID/build"
mkdir -p "$builddir"
cd "$builddir"
CC=gcc-7 CXX=g++-7 cmake ~/gromacs-constantph -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=${PWD}/.. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA
make -j 12
make install -j 12
cd ..
source ${PWD}/bin/GMXRC

# RUN SIMULATION

cd "$simdir"

if [ ! -e MD.tpr ]
then
    gmx grompp -f MD.mdp -c CA.gro -p topol.top -n index.ndx -o MD.tpr -maxwarn 1
    gmx mdrun -deffnm MD -c MD.pdb -x MD.xtc -npme 0 -maxh 47 -bonded gpu -nt $SLURM_JOB_CPUS_PER_NODE
else
    gmx mdrun -deffnm MD -c MD.pdb -x MD.xtc -cpi MD.cpt -npme 0 -maxh 47 -bonded gpu -nt $SLURM_JOB_CPUS_PER_NODE
fi

# RESUBMIT IF NOT DONE

if [ ! -e MD.pdb ]
then
    sbatch jobscript.sh
fi
