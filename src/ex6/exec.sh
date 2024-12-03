#bin/bash
set -e
gcc -o ex6 ex6.c MathTool.c Mesh.c SolidFEM.c GLTool.c -lglut -lGLU -lGL -lm
./ex6