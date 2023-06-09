#! /bin/bash

i=131
eqS=5
S=40
N=3

sysname=`cat name`
NREPS=`cat nreps`

DIR=`pwd`
PANADIR=$DIR/analysis$(( $i - 1 ))
ANADIR=$DIR/analysis${i}

if [ ! -d $ANADIR ]; then
  mkdir $ANADIR
fi
if [ ! -d $ANADIR/data ]; then
  mkdir $ANADIR/data
fi

cd $ANADIR

../ALF/GetVolumes.py $i $N $eqS $S
