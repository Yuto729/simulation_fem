#bin/bash
set -e
gcc -o ex4 ex4.c MathTool.c Mesh.c SolidFEM.c GLTool.c -lglut -lGLU -lGL -lm
./ex4