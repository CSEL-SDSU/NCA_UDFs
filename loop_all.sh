#!/bin/bash

# To run: 
#   ./loop_all.sh [base case file] 
#   ./loop_all.sh [base case file] my_suffix
# inside the directory with the case file and the NCA_unix.c file
# make executable by doing chmod +x submit_fluent.sh

# This script is used to queue jobs in Earth, Lunar, and Martial gravity
# for horizontal and vertial configurations

# Rule for inputs ${var:-default} 
BASE_CASE=${1:?Usage: ./loop_all.sh base_case.cas.h5 [append]}
APPEND=${2:-}
H=${3:-0.00495} #default gap height of 0.00495m

if [[ -n "$APPEND" ]]; then
    APPEND_SUFFIX="_$APPEND"
else
    APPEND_SUFFIX=""
fi

# Add jobs for horizontal and vertical burns in Earth, Lunar, and Martian gravity
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

        GX="$GX" GY="$GY" H="$H" ./submit_fluent_2.sh "$RUN_NAME" "$BASE_CASE"
        #echo "$RUN_NAME"
    done
done

# Add jobs for adittional angles -75:15:75 degrees in Earth gravity
ANGLES=($(for ((i=75; i>=-75; i-=15)); do echo "$i"; done))

#Loop over angles
for angle in "${ANGLES[@]}"; do
    GX=$(awk -v ang="$angle" -v G="$G" -v PI="$PI" 'BEGIN {print G * sin(ang * PI / 180.0)}')
    GY=$(awk -v ang="$angle" -v G="$G" -v PI="$PI" 'BEGIN {print G * cos(ang * PI / 180.0)}')

    RUN_NAME="tilt_${angle}deg_${APPEND_SUFFIX}"

    GX="$GX" GY="$GY" H="$H" ./submit_fluent_2.sh "$RUN_NAME" "$BASE_CASE"
    #echo "GX=""$GX"
    #echo "GY=""$GY"
done

