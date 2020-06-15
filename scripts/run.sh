#!/bin/bash
# Run all X jobs on Y CPUs (2*Y vCPUs)
# Do not forget to time it
# parallel -j Y ./task_runner_jsp.sh {%} {} ::: {1..X}

parallel -j 1 ./task_runner_jsp.sh {%} {} ::: {1..1}
