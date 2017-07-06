#!/bin/bash
# Script for quick compilation of compuGUT source file
#
# Programmer: Arun S. Moorthy, amoorthy@uoguelph (canada)
# Revision: 0.3
# Date: August 19th, 2015
# ------------------------------------------------------------------------------

# compile from source (using gcc)
gcc compuGUTv01.c -I ./SUNDIALS_LIB/include/ \
                  -L ./SUNDIALS_LIB/lib/ \
                  -lsundials_kinsol\
                  -lsundials_cvode \
                  -lsundials_nvecserial \
                  -lm \
                  -o compuGUTv01-l64g

# clean folder of temp files
rm *.*~

# clear screen
# clear
