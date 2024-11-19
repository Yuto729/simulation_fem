#!/bin/bash
set -e

gcc -o ex1 ex1.c MathTool.c -lm
./ex1 > result.dat
diff -u result.dat answer.dat
