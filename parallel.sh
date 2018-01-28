#!/bin/bash
#SBATCH -N 1  # 1 node
#SBATCH -n 1  # 1 task per node
#SBATCH -c 16 # Number of cores per task.
#SBATCH -o matlab.out # stdout is redirected to that file
#SBATCH -e matlab.err # stderr is redirected to that file
#SBATCH --gres=gpu:1 # Use GPUs
#SBATCH -t 20:00:00 
# sends mail when process begins, and
# when it ends. Make sure you define your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
# replace your_netid by your actual netid
#SBATCH --mail-user=eboix@princeton.edu

srun /usr/licensed/bin/matlab -nodesktop -nosplash -r "addpath('/home/eboix/sum17');conn_stats;quit"
