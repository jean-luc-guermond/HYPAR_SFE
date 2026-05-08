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
    mv mesh_name previous_data* data_regression* mesh_part* Mesh_1* *.plt job.sh output_$i

    #===Clean up
    rm -rf output_$i
done
