#!/bin/bash
# Run all N jobs on C CPUs (2*C vCPUs)
# Do not forget to time it
# parallel -j C ./task_runner.sh {%} {} ::: {1..N}

parallel -j 1 ./task_runner.sh {%} {} ::: {1..1}
