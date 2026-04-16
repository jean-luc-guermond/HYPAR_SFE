#!/bin/bash
datatest=data

nproc=2

for ((i=1; i<=$3; i++)); do
    #=== define executable
    exe_index=$((3 + i))
    exe=${!exe_index}
    #=== run the test
    $1 $2$nproc ../EXECUTABLE/${exe} regression $i
    #=== move the output
    mkdir output_$i
    mv previous_data mesh_part* Mesh_1* *.plt output_$i
done

#Clean up
rm -rf output*
