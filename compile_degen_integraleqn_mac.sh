#!/bin/bash

g++-13 ./source/efimov_code_check.cpp -O3 -std=c++14 -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/ -L/usr/local/lib -fopenmp -o efimov_printer_macos
