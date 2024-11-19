#!/bin/bash
set -e

gcc -o ex2 ex2.c MathTool.c Mesh.c SolidFEM.c -lm
./ex2>result.dat
diff -u -w -B result.dat answer.dat