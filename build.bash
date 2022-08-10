#!/bin/bash

mkdir -p bin

for i in {1..12}
do
    echo "Compiling step $i..."
    chpl --fast -o bin/step_$i chapel_ports/step_$i.chpl
done

for i in {4..12}
do
    echo "Compiling step $i distributed..."
    chpl --fast -o bin/step_${i}_dist chapel_ports/step_${i}_dist.chpl
done
