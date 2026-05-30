#!/bin/bash
i_nproc=$((3 + $3 + 1))
nproc=${!i_nproc}

for ((i=1; i<=$3; i++)); do
    #=== define executable
    exe_index=$((3 + i))
    exe=${!exe_index}
    
    #=== run the test
    $1 $2$nproc ../EXECUTABLE/${exe} regression $i

    #=== move the output
    mkdir output_$i
    shopt -s nullglob #to avoid warning msg
    mv previous_data* data_regression* mesh_name mesh_part* Mesh_1* *.plt output_$i
    shopt -u nullglob #to avoid warning msg

    #===Clean up
    rm -rf output_$i
done
