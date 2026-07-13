#!/bin/bash
i_nproc=$((3 + $3 + 1))
nproc=${!i_nproc}

for ((i=1; i<=$3; i++)); do
    #=== define executable
    exe_index=$((3 + i))
    exe=${!exe_index}
    
    #=== run the test
    $1 $2$nproc ../EXECUTABLE/${exe} if_regression=.True. num_regex=${i} \
    data_in=data_${i} \
    data_save=previous_data_${i} \
    data_out=data_regression_${i}_NPROC_${nproc}
    #=== run the test

    #=== run the test while rewriting data files
    ##=== UNCOMMENT WITH PRECAUTION ===#
    # $1 $2$nproc ../EXECUTABLE/${exe} if_regression=.True. num_regex=${i} \
    # data_in=data_${i} \
    # data_save=previous_data_${i} \
    # data_out=data_${i}
    ##=== UNCOMMENT WITH PRECAUTION ===#

    #=== move the output
    mkdir output_$i
    shopt -s nullglob #to avoid warning msg
    mv previous_data* data_regression* mesh_name mesh_part* Mesh_1* *.plt backup* series* *.pvd *.vtu output_$i 
    shopt -u nullglob #to avoid warning msg

    #===Clean up the output
    rm -rf output_$i
done
