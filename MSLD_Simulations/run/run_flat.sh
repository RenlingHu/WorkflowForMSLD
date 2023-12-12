#!/bin/sh

cd ../

export nnodes=`cat nnodes`
export nreps=`cat nreps`

# DEPEND="--dependency=afterok:"
bash runset2.sh
PID=`sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset2.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
PID=`sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset3.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
PID=`sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset3_2.sh | awk '{print $4}'`
echo $PID

done