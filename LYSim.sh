#!/bin/bash

#
# variables from arguments string in jdl
#

RUN_DIR=$1
INPUT_NAME=$2
ANALYSIS_FILE_NAME=$3

#
# header 
#
cd $RUN_DIR
source /data/users/yhshin/geant4/geant4.9.6.p02-install/bin/geant4.sh

./LYSim ${INPUT_NAME} ${ANALYSIS_FILE_NAME}
