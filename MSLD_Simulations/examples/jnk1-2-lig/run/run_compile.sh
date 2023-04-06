#!/bin/sh

module load cuda/11.2
module load gcc

cd ../dWHAMdV_mso/
bash Clean.sh

bash Compile.sh