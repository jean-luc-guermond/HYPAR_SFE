#!/bin/bash
datatest=data                                    
                                                             
nproc=32

$1 $2$nproc $4  ../EXECUTABLE/$3 regression
echo $?
#cp current_regression_reference regression_reference_10

#Clean up
rm -f previous_data mesh_part* Mesh_1* *.plt
