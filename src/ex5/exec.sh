#bin/bash
set -e
gcc -o ex5 ex5.c MathTool.c Mesh.c SolidFEM.c GLTool.c -lglut -lGLU -lGL -lm
./ex5