#! /bin/bash

cd ../

export nnodes=`cat nnodes`
export nreps=`cat nreps`

# phase 1-3
#DEPEND="--dependency=afterok:$PID"
PID=`sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset2.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
PID=`sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset3.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
PID=`sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset3_2.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"


# production run
export nitt=10

NEWDEPEND="--dependency="
COMMA=""
for p in a b c
do

export ini=131
export i=$ini$p
eval PID=`sbatch --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:$nnodes --export=ALL --array=1-4%1 $DEPEND ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","
echo "$p-$PID"

done

# Data analysis
DEPEND=$NEWDEPEND
export i=131
export eqS=5
export S=40
export N=3
# export skipE=10

PID=`sbatch -N 1 --ntasks=$nreps --cpus-per-task=$nnodes -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"