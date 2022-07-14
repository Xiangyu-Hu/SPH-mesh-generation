#!/bin/sh

filename=output.msg

echo MRSPH code is started! | tee $filename
echo Start time is `date` | tee -a $filename

make -j run | tee -a $filename
# make -j

# mpiexec -n 1 -host localhost ./MRSPH | tee -a $filename

echo End time is `date` | tee -a $filename

mv $filename ./outdata/

