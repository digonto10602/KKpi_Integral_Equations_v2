#!/bin/bash

g++ ./source/efimov_code_check.cpp -O3 -std=c++14 -I/usr/include/eigen3/ -fopenmp -o efimov_printer
