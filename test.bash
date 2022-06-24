#!/bin/bash

for i in {1..4}
do
    echo "Comparing Chapel and Python results from step $i"
    chpl --fast -o bin/step_$i chapel_ports/step_$i.chpl
    bin/step_$i > /dev/null
    python3 python_scripts/step_$i.py > /dev/null
    diff "sim_output/step_${i}_output.txt" "sim_output/step_${i}_py_output.txt"
done
