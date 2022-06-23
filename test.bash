#!/bin/bash

for i in {1..4}
do
    echo "Comparing Chapel and Python results from step $i"
    chpl ./step_$i.chpl
    ./step_$i > /dev/null
    python3 step_$i.py > /dev/null
    diff "sim_output/step_${i}_output.txt" "sim_output/step_${i}_py_output.txt"
done
