#!/bin/bash
datatest=data                                    
                                                             
nproc=32

$1 $2$nproc $4  ../EXECUTABLE/$3 regression
echo $?

#Clean up
rm -f *.plt previous_data
