#!/bin/bash

# To run: 
#   ./loop_grav_accel.sh [base case file] 
#   ./loop_grav_accel.sh [base case file] my_suffix
# inside the directory with the case file and the NCA_unix.c file
# make executable by doing chmod +x submit_fluent.sh

# This script is used to queue jobs in Earth, Lunar, and Martial gravity
# for horizontal and vertial configurations

BASE_CASE=${1:?Usage: ./loop_grav_accel.sh base_case.cas.h5 [append]}
APPEND=${2:-}

if [[ -n "$APPEND" ]]; then
    APPEND_SUFFIX="_$APPEND"
else
    APPEND_SUFFIX=""
fi

ANGLES=("0" "90" "-90")
G_BY_GE=("1" "0.166" "0.38")

GE=-9.81
PI=3.1415927

for grav_frac in "${G_BY_GE[@]}"; do
    for angle in "${ANGLES[@]}"; do

        G=$(awk -v g_ov_ge="$grav_frac" -v ge="$GE" 'BEGIN {print g_ov_ge * ge}')

        GX=$(awk -v ang="$angle" -v G="$G" -v PI="$PI" 'BEGIN {print G * sin(ang * PI / 180.0)}')
        GY=$(awk -v ang="$angle" -v G="$G" -v PI="$PI" 'BEGIN {print G * cos(ang * PI / 180.0)}')

        
        RUN_NAME="tilt_${angle}deg_gfrac_${grav_frac}${APPEND_SUFFIX}"

        GX="$GX" GY="$GY" ./submit_fluent_grav.sh "$RUN_NAME" "$BASE_CASE"
        #echo "$RUN_NAME"
    done
done


