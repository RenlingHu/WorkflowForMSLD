#! /bin/bash

cd ../

export nnodes=`cat nnodes`
export nreps=`cat nreps`

export i=131
export eqS=5
export S=40
export N=3

sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL ./postvolume.sh

done