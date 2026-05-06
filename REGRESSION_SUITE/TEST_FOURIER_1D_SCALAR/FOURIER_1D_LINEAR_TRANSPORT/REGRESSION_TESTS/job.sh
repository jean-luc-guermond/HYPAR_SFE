#!/bin/bash
datatest=data

# nproc=1
i_nproc=$((3 + $3 + 1))
nproc=${!i_nproc}

for ((i=1; i<=$3; i++)); do
    #=== define executable
    exe_index=$((3 + i))
    exe=${!exe_index}
    #=== run the test
    cp data_$i data
    $1 $2$nproc ../EXECUTABLE/${exe} regression $i
    #=== move the output
    mkdir -p output_$i
    mv previous_data *.plt output_$i
done

#Clean up
rm -f data
# rm -rf output_*
