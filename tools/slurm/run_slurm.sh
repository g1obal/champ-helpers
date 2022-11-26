#!/bin/bash

# run_slurm
#
# Copy slurm files to run directories.
# Push to run directory, sbatch slurm jobs.
#
# Author: Gokhan Oztarhan
# Created date: 16/02/2022
# Last modified: 01/11/2022

ROOT_DIR="."
INPUT_FILE_NAME="input_file"

# Find directories including INPUT_FILE_NAME recursively.
# Does not work if directory names include white spaces or special characters!
paths=()
for D in `find $ROOT_DIR -type f -name $INPUT_FILE_NAME | sort -V`; do
    paths+=(${D%/*}) # run paths, removing INPUT_FILE_NAME from string
done

# Loop over paths
for i in ${!paths[@]}; do
    run_path=${paths[$i]}
    
    cp "slurm_$HOSTNAME" $run_path
    
    pushd $run_path > /dev/null
    
    sbatch "slurm_$HOSTNAME" > run_info
    echo $run_path
    
    popd > /dev/null
done

