#! /bin/bash
# task_runner_list.sh
 
slot=$1
 
jobfile=$2
 
# This function determines which virtual CPUs correspond to the specified slot:
cpu_threads() {
    minus=$((slot - 1))
    first_cpu_thread=$((2 * minus))
    second_cpu_thread=$((first_cpu_thread + 1))
    echo $first_cpu_thread,$second_cpu_thread
}

taskset -c $(cpu_threads) python runTractor.py ${jobfile}
