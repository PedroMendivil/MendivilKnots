#!/bin/bash

# Ejecuta el programa simulate-polymer.bin con el directorio como argumento

#                               ./simulate-polymer.bin Simulations/Test4

# Loop para ejecutar vmd con diferentes archivos de trayectoria
for i in {0..2}; do
    vmd -e visualize-polymer.tcl -args Simulations/Test4/initial-condition.gro Simulations/Test4/trajectory-file-$i.trr
done

