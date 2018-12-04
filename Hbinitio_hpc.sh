#! /bin/bash
####SBATCH -p public
#SBATCH -p grant -A g2018a7
###SBATCH -p grantgpu -Ag1018a7
#SBATCH -n 64 -tasks-per-node=64
#SBATCH -t 2:00:00                # Le job sera tue au bout de 1h
####SBATCH --mail-user=herve.bulou@ipcms.unistra.fr
#####SBATCH --mem=4096                # Quantite memoire demandee par noeud en Mo (unite obligatoire)


export OMP_NUM_THREADS=64
time ./Hbinitio.x

