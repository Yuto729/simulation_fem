#!/bin/bash

# コンパイル
gcc -o ex6 ex6.c SolidFEM.c GLTool.c MathTool.c Mesh.c -lglut -lGLU -lGL -lm

# 実行
if [ $? -eq 0 ]; then
    ./ex6
else
    echo "Compilation failed"
fi 