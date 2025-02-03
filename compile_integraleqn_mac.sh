#!/bin/bash

#rm printer.o

g++-13 ./source/printer_test.cpp -O3 -std=c++14 -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/ -fopenmp -o printer_macos
