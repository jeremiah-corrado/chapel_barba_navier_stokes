#!/bin/bash

for i in {1..9}
do
    echo "Comparing Chapel and Python results from step $i..."
    mkdir -p sim_output/step_$i
    chpl --fast -o bin/step_$i chapel_ports/step_$i.chpl
    bin/step_$i > /dev/null
    python3 python_scripts/step_$i.py > /dev/null
    diff "sim_output/step_${i}/ch_u.txt" "sim_output/step_${i}/py_u.txt"
done

for i in {4..9}
do
    echo "Comparing Chapel Distributed and Python results from step $i..."
    chpl --fast -o bin/step_$i chapel_ports/step_${i}_distributed.chpl
    bin/step_$i --write_data=true > /dev/null
    diff "sim_output/step_${i}/ch_u.txt" "sim_output/step_${i}/py_u.txt"
done
