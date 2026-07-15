#!/usr/bin/env bash

# Clear any conflicting modules
module purge

# Load the GCC 13.3.0 toolchain and dependencies for cneta
module load GCC/13.3.0
module load CMake/3.29.3
module load Boost/1.85.0 
module load GSL/2.8

echo "HPC environment loaded successfully for cneta."
echo "You can now run: cd ../code/build && cmake .. && make -j4"
