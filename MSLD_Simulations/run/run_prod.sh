#! /bin/bash


cd ../

export nnodes=`cat nnodes`
export nreps=`cat nreps`
export nitt=10

for p in a b c
do

export ini=131
export i=${ini}$p
sbatch --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:$nnodes --export=ALL --array=1-4%1 ./runset4.sh

done