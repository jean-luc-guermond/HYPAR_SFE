#!/bin/bash
datatest=data

nproc=4

for ((i=1; i<=$3; i++)); do
    exe_index=$((3 + i))
    exe=${!exe_index}
echo "$1 $2$nproc  ../EXECUTABLE/${exe} regression $i"
    $1 $2$nproc ../EXECUTABLE/${exe} regression $i
done
#cp current_regression_reference regression_reference_10

#Clean up
rm -f previous_data mesh_part* Mesh_1* *.plt
