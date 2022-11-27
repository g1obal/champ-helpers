"""
job_queue

Run CHAMP jobs using mpirun with a queue system.

BEFORE USE:
Modify CHAMP variable for your CHAMP path.

If QUEUE=True and WAIT=True,
    This script finds all directories containing the input file
    named by INPUT_FILE_NAME variable, push directory and runs CHAMP.
    Waits for job to finish before executing another. Runs single job at a time.

If QUEUE=True and WAIT=False,
    Does not wait for jobs. Executes all runs simultaneously.

CHAMP run commands:
mpirun -np 2 champ_mpi.exe -m vmc_mov1_mpi -i input >output
mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi1 -i input >output
mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi2 -i input >output
mpirun -np 2 champ_mpi.exe -m dmc_mov1_mpi3 -i input >output

Author: Gokhan Oztarhan
Created date: 15/01/2020
Last modified: 27/11/2022
"""

import os
import errno
import time
from datetime import datetime
from subprocess import Popen, PIPE


N_CPU = 8
RUN_MODE = 'vmc'

QUEUE = True
ROOT_DIR = '.'
WAIT = True
SLEEP = 5

INPUT_FILE_NAME = 'input_file'
OUTPUT_FILE_NAME = 'output_file'

MPIRUN = '/usr/bin/mpirun' # use '$ which mpirun' to find mpirun path
HOME = os.path.expanduser('~')
CHAMP = os.path.join(HOME, 'programs/CHAMP/qmc/champ_mpi.exe')
MODE = {
    'vmc': '-m vmc_mov1_mpi',
    'dmc1': '-m dmc_mov1_mpi1',
    'dmc2': '-m dmc_mov1_mpi2',
    'dmc3': '-m dmc_mov1_mpi3',
}


def run_champ():
    run_command = ' '.join([
        MPIRUN, '-n %d --use-hwthread-cpus --bind-to none' %(N_CPU), 
        CHAMP, MODE[RUN_MODE],
        '-i', INPUT_FILE_NAME, 
        '>', OUTPUT_FILE_NAME, 
        '&', 'echo $!',
    ])
    
    sub = Popen(run_command, shell=True, stdout=PIPE)
    PID_sub = sub.pid # sub shell pid
    PID_champ = sub.stdout.readline().decode() # echo $!
    sub.stdout.close()
    sub.kill()

    with open('run_info', 'w') as f:
        f.write(
            'run_command = %s\n\n' %run_command \
            + 'PID_python = %s\n' %os.getpid() \
            + 'PID_sub = %s\n' %PID_sub \
            + 'PID_champ = %s\n' %PID_champ
        )
    
    return int(PID_champ)


def job_queue():
    # Find directories including INPUT_FILE_NAME recursively.
    paths = [
        root for root, dirs, files in sorted(os.walk(ROOT_DIR)) \
        if INPUT_FILE_NAME in files
    ]
    
    # Log file
    log_file = open('log-job_queue','w')
    
    # Print info
    start_time = datetime.now()
    log_file.write(
        'PID_python = %s\n' %os.getpid() \
        + 'start_time = %s\n\n' %start_time.strftime('%Y-%m-%d %H:%M:%S') \
        + 'N_CPU = %i\n' %N_CPU \
        + 'RUN_MODE = %s\n\n' %RUN_MODE \
        + 'QUEUE = %s\n' %QUEUE \
        + 'ROOT_DIR = %s\n' %ROOT_DIR \
        + 'WAIT = %s\n' %WAIT \
        + 'SLEEP = %s\n\n' %SLEEP \
        + 'total_runs = %s\n\n' %len(paths) \
        + '%-7s %-7s %-s\n' %('index', 'PID', 'run_path') \
        + '%-7s %-7s %-s\n' %('-----', '---', '--------')
    )
    
    # Working directory
    root = os.getcwd()
    
    for i, run_path in enumerate(paths):
        # push to run_directory
        os.chdir(run_path)
        
        # Run CHAMP
        PID_champ = run_champ()
        
        # Print info
        log_file.write('%-7s %-7s %-s\n' %(i, PID_champ, run_path))
        log_file.flush()
        
        # Wait for finish if desired
        if WAIT:
            while is_running(PID_champ):
                time.sleep(SLEEP)
        
        # pop back to working directory
        os.chdir(root)
        
    # Print info
    end_time = datetime.now()
    log_file.write(
        '\nend_time = %s\n' %end_time.strftime('%Y-%m-%d %H:%M:%S') \
        + 'elapsed_time = %s\n' %str(end_time - start_time)
    )
    
    # Close log file
    log_file.close()

    
def is_running(pid):     
    """
    Check if pid is running using os.kill(pid, 0).
    
    Bash equivalent: 'kill -0 pid' 
    '-0' flag just checks whether pid is running and accesible or not.
    It does not send any signal.
    
    No Error: The process exists.
    errno.EPERM: Operation not permitted.
                 The process exists but does not belong to user.
    errno.ESRCH: No such process.
    
    Source: https://stackoverflow.com/a/7654102/13893858
    """   
    try:
        os.kill(pid, 0)
    except OSError as err:
        if err.errno == errno.ESRCH:
            return False
    return True
    

if __name__ == '__main__':
    if QUEUE:
        job_queue()
    else:
        run_champ()


