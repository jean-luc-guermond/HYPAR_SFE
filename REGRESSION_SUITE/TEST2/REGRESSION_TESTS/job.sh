#!/bin/bash
datatest=data                                    
                                                             
nproc=4

$1 $2$nproc $4  ../EXECUTABLE/$3 regression
echo $?
#cp current_regression_reference regression_reference_10

#Clean up
rm -f *.plt previous_data