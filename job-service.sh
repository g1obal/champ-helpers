#!/bin/bash

# job-service
#
# Run jobs using simple queue system.
#
# BEFORE USE MODIFY:
# N_CORES_PER_JOB, RUN_MODE, 
# INPUT_FILE_NAME, OUTPUT_FILE_NAME, ERROR_FILE_NAME, CHAMP
#
# FOR INPUT_FILE SEARCH:
# Directory names should not include spaces or special characters!
#
# jobs_waiting: The file containing waiting jobs. Jobs can be added or removed
#               while the service is running.
# jobs_running: The file containing currently running jobs.
# jobs_log: Log file including finished jobs and service information.
#
# Author: Gokhan Oztarhan
# Created date: 27/05/2024
# Last modified: 29/05/2024

N_CORES_PER_JOB=2

RUN_MODE="vmc_mov1_mpi" # vmc_mov1_mpi, dmc_mov1_mpi1

INPUT_FILE_NAME="input_file"
OUTPUT_FILE_NAME="output_file"
ERROR_FILE_NAME="error_file"

CHAMP="$HOME/programs/CHAMP/qmc/champ_mpi.exe"

SLEEP_INTERVAL=2

# Determine the number of physical cores
MAX_CORES=$(lscpu | awk '
    /^Core\(s\) per socket:/ {cores_per_socket=$4}
    /^Socket\(s\):/ {sockets=$2}
    END {print cores_per_socket * sockets}
')

# Determine the max number of jobs
MAX_JOBS=$(( MAX_CORES / N_CORES_PER_JOB ))

# core: bind processes to cores, none: OS scheduler distributing processes
MPI_BIND_TO="none"

# Terminate the script, if CHAMP file does not exist.
if [ ! -f $CHAMP ]; then
    echo "CHAMP is not found."
    exit
fi

submit () {    
    # Open empty file (overrides existing file)
    > jobs_waiting

    # Find directories including INPUT_FILE_NAME recursively from pwd.
    # Save run real path (remove INPUT_FILE_NAME from the string), and command
    for D in `find . -type f -name $INPUT_FILE_NAME | sort -V`; do
        echo "`realpath ${D%/*}`," \
            "mpirun --bind-to $MPI_BIND_TO -np $N_CORES_PER_JOB" \
            "$CHAMP -m $RUN_MODE -i $INPUT_FILE_NAME" \
            "1> $OUTPUT_FILE_NAME 2> $ERROR_FILE_NAME" >> jobs_waiting
    done
}

service () {
    # Open empty file (overrides existing file)
    > jobs_log
    
    # Remove if there is jobs_running file
    if [ -f jobs_running ]; then
        rm jobs_running
    fi
    > jobs_running
    
    # Print info
    printf "PID service = %s\n" $$ >> jobs_log
    printf "start time = %s %s\n\n" `date +"%F %T"` >> jobs_log
    printf "N_CORES_PER_JOB = %s\n" $N_CORES_PER_JOB >> jobs_log
    printf "MAX_CORES = %s\n" $MAX_CORES >> jobs_log
    printf "MAX_JOBS = %s\n\n" $MAX_JOBS >> jobs_log
    printf "JOBS FINISHED (PID, path, command)\n" >> jobs_log
    echo "----------------------------------" >> jobs_log
    
    # PID array of running jobs
    declare pid_running
    
    n_jobs_running=0
    while true; do
        # If there no lines in jobs_waiting, break the loop
        if [ ! -s jobs_waiting ]; then 
            break
        fi
        
        # Wait for running jobs
        while [ $n_jobs_running -ge $MAX_JOBS ]; do
            for pid in "${!pid_running[@]}"; do
                if ! kill -0 $pid 2> /dev/null; then
                    # Job finished
                    (( n_jobs_running -= 1 ))
                    unset pid_running["$pid"]
                    
                    # Get entry from jobs_running, add it to jobs_log
                    echo $(grep "^$pid," jobs_running) >> jobs_log
                    
                    # Remove the entry from jobs_running
                    sed -i "/^$pid,/d" jobs_running
                fi
            done
            
            sleep $SLEEP_INTERVAL
        done
        
        # Read path and command from jobs_waiting in the 1st line
        IFS="," read -r path command < jobs_waiting
        
        # Remove the first line from jobs_waiting
        sed -i "1d" jobs_waiting
        
        # Create a temporary named pipe for PID storage
        pid_temp=$(mktemp -u /tmp/temp_pid.XXXXXX)
        mkfifo $pid_temp
        
        # Run the command in a new shell
        bash -c "cd $path; $command & echo \$! > $pid_temp & disown -a; exit" &
        
        # Read PID and remove temporary pipe
        read pid < $pid_temp
        rm $pid_temp
        
        # Add entry to jobs_running
        echo "$pid, $path, $command" >> jobs_running
        
        # Add element to pid_running array
        pid_running["$pid"]=$pid
        
        # Increase job running count
        (( n_jobs_running ++ ))
        
        sleep $SLEEP_INTERVAL
    done
    
    # Wait for all runs to finish
    wait
    
    # Remaining jobs
    while [ $n_jobs_running -ge 1 ]; do
        for pid in "${!pid_running[@]}"; do
            if ! kill -0 $pid 2> /dev/null; then
                # Job finished
                (( n_jobs_running -= 1 ))
                unset pid_running["$pid"]
                
                # Get entry from jobs_running, add it to jobs_log
                echo $(grep "^$pid," jobs_running) >> jobs_log
                
                # Remove the entry from jobs_running
                sed -i "/^$pid,/d" jobs_running
            fi
        done
        
        sleep $SLEEP_INTERVAL
    done
    
    # Remove service files
    rm jobs_waiting
    rm jobs_running

    # Print info
    printf "\nend time = %s %s\n" `date +"%F %T"` >> jobs_log
    secs=$SECONDS
    hrs=$(( secs/3600 ))
    mins=$(( (secs-hrs*3600)/60 ))
    secs=$(( secs-hrs*3600-mins*60 ))
    printf "elapsed_time = %02d:%02d:%02d\n" $hrs $mins $secs >> jobs_log
}

waiting () {
    echo $(wc -l jobs_waiting)
}

running () {
    echo $(wc -l jobs_running)
}

stop () {
    kill $(awk '/PID service =/ {pid_service=$4} 
        END {print pid_service}' jobs_log)
    IFS=","
    while read -r pid path command
    do
        kill $pid
    done < jobs_running
}

export N_CORES_PER_JOB RUN_MODE
export INPUT_FILE_NAME OUTPUT_FILE_NAME ERROR_FILE_NAME CHAMP
export SLEEP_INTERVAL MAX_CORES MAX_JOBS MPI_BIND_TO
export -f submit
export -f service

case "$1" in
    submit)
        bash -c 'submit' &
        ;;
    service)
        bash -c 'service' &
        ;;
    ss)
        submit &&
        bash -c 'service' &
        ;;
    waiting)
        waiting
        ;;
    running)
        running
        ;;
    stop)
        stop
        ;;
    *)
        echo "usage $0 {submit,service,ss,waiting,running,stop}"
        exit 1
        ;;
esac

exit 0


