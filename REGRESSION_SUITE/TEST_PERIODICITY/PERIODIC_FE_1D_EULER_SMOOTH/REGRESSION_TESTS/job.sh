#!/bin/bash
datatest=data

nproc=8

for ((i=1; i<=$3; i++)); do
    #=== define executable
    exe_index=$((3 + i))
    exe=${!exe_index}
    #=== run the test
    $1 $2$nproc ../EXECUTABLE/${exe} regression $i
    #=== move the output
    mkdir output_$i
    mv previous_data *.plt output_$i
done

#Clean up
rm -rf output*