# bin/bash
set -e
gcc test.c -o test -lglut -lGLU -lGL
./test