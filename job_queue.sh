#!/bin/bash

# job_queue
#
# Run CHAMP jobs using simple queue system.
#
# BEFORE USE:
# Modify CHAMP variable for your CHAMP path.
# 
# MAX_JOBS is the maximum number of jobs to run simultaneously.
# MAX_JOBS * N_CPU_PER_JOB should be lower than (or equal to) the maximum number
# of cores available in the system !!
#
# Directory names should not include spaces or special characters !!
#
# CHAMP run commands:
# mpirun -np 2 champ_mpi.exe -m vmc_mov1_mpi -i input > output
# mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi1 -i input > output
# mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi2 -i input > output
# mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi3 -i input > output
#
# Author: Gokhan Oztarhan
# Created date: 19/02/2022
# Last modified: 18/01/2023

N_CPU_PER_JOB=8
MAX_JOBS=5

RUN_MODE="vmc"

ROOT_DIR="."
SLEEP_INTERVAL=5

INPUT_FILE_NAME="input_file"
OUTPUT_FILE_NAME="output_file"
ERROR_FILE_NAME="error_file"

CHAMP="$HOME/programs/CHAMP/qmc/champ_mpi.exe"

declare -A MODE=(
    ["vmc"]="-m vmc_mov1_mpi"
    ["dmc1"]="-m dmc_mov1_mpi1"
    ["dmc2"]="-m dmc_mov1_mpi2"
    ["dmc3"]="-m dmc_mov1_mpi3"
)

# Terminate the script, if CHAMP file does not exist.
if [ ! -f $CHAMP ]; then
    echo "CHAMP is not found."
    exit
fi

# Find directories including INPUT_FILE_NAME recursively.
# Does not work if directory names include white spaces or special characters!
paths=()
for D in `find $ROOT_DIR -type f -name $INPUT_FILE_NAME | sort -V`; do
    paths+=(${D%/*}) # run paths, removing INPUT_FILE_NAME from string
done

# Open file and assign 10 as file handler
exec 10> log-job_queue

# Print info
printf "PID_bash = %s\n" $$ >& 10
printf "start_time = %s %s\n\n" `date +"%F %T"` >& 10
printf "N_CPU_PER_JOB = %s\n" $N_CPU_PER_JOB >& 10
printf "MAX_JOBS = %s\n\n" $MAX_JOBS >& 10
printf "RUN_MODE = %s\n\n" $RUN_MODE >& 10
printf "ROOT_DIR = %s\n" $ROOT_DIR >& 10
printf "SLEEP_INTERVAL = %s\n\n" $SLEEP_INTERVAL >& 10
printf "total_runs = %s\n\n" ${#paths[@]} >& 10
printf "%-7s %-7s %-s\n" "index" "PID" "run_path" >& 10
printf "%-7s %-7s %-s\n" "-----" "---" "--------" >& 10

# Main loop for jobs
ipaths=0
while (( $ipaths < ${#paths[@]} )) ; do
    # Run new job if number of running jobs is lower than MAX_JOBS
    if (( `jobs -rp | wc -l` < $MAX_JOBS )) ; then
        pushd ${paths[$ipaths]} > /dev/null

        # Open MPI "--bind-to core" is the default setting. If more than one MPI
        # jobs are started from the same bash script or terminal, the processes
        # of all MPI jobs (indicated by -np flag) share the first N cores. Thus, 
        # time sharing of the cores slows the execution of MPI processes.
        # Use "--bind-to none" to let the Linux scheduler use all the available
        # cores. Source: https://stackoverflow.com/a/66112173/13893858
        mpirun -np $N_CPU_PER_JOB --use-hwthread-cpus --bind-to none \
        $CHAMP ${MODE[$RUN_MODE]} -i $INPUT_FILE_NAME 1> $OUTPUT_FILE_NAME \
        2> $ERROR_FILE_NAME &

        PID_champ=$!
        echo "PID_bash = $$" > run_info
        echo "PID_champ = $PID_champ" >> run_info

        printf "%-7s %-7s %-s\n" $ipaths $PID_champ ${paths[$ipaths]} >& 10
        popd > /dev/null
        
        ipaths=$(( ipaths + 1 ))
    fi
    
    sleep $SLEEP_INTERVAL
done

# Wait for all runs to finish
wait

# Print info
printf "\nend time = %s %s\n" `date +"%F %T"` >& 10
secs=$SECONDS
hrs=$(( secs/3600 ))
mins=$(( (secs-hrs*3600)/60 ))
secs=$(( secs-hrs*3600-mins*60 ))
printf "elapsed_time = %02d:%02d:%02d\n" $hrs $mins $secs >& 10

# Close file
exec 10<&-


