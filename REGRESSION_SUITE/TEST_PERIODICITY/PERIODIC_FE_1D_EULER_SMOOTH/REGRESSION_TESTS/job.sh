#!/bin/bash
datatest=data

nproc=8
$1 $2$nproc $4  ../EXECUTABLE/$3 regression
#cp current_regression_reference regression_reference_10

#Clean up
#rm -f previous_data *.plt
