#!/bin/bash

# To run: ./loop_grav.sh [base case file] 
# inside the directory with the case file and the NCA_unix.c file
# make executable by doing chmod +x submit_fluent.sh

BASE_CASE=${1:?Usage: ./submit_fluent.sh base_case.cas.h5}

#Generate angles
# For all angles:
ANGLES=($(for ((i=60; i>=-75; i-=15)); do echo "$i"; done))

#For one angle (debugging):
#ANGLES=75

G=-9.81
PI=3.1415927

#Loop over angles
for angle in "${ANGLES[@]}"; do
    GX=$(awk -v ang="$angle" -v G="$G" -v PI="$PI" 'BEGIN {print G * sin(ang * PI / 180.0)}')
    GY=$(awk -v ang="$angle" -v G="$G" -v PI="$PI" 'BEGIN {print G * cos(ang * PI / 180.0)}')

    RUN_NAME=tilt_"$angle"deg

    GX="$GX" GY="$GY" ./submit_fluent_grav.sh "$RUN_NAME" "$BASE_CASE"
    #echo "GX=""$GX"
    #echo "GY=""$GY"
done