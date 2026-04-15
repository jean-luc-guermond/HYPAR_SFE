#!/bin/bash
datatest=data                                    
                                                             
nproc=1

$1 $2$nproc $4  ../EXECUTABLE/$3 regression
echo $?

#Clean up
rm -f previous_data *.plt