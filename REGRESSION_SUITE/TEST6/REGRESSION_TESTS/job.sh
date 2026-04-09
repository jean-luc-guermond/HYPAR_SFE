#!/bin/bash
datatest=data                                    
                                                             
nproc=1

$1 $2$nproc $4  ../EXECUTABLE/$3
echo $?

#Clean up
rm -f previous_data mesh_part* Mesh_1* *.plt
