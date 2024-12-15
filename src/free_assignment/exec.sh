#!/bin/bash

# コンパイル
gcc -o main main.c SolidFEM.c GLTool.c MathTool.c Mesh.c -lglut -lGLU -lGL -lm

# 実行
if [ $? -eq 0 ]; then
    ./main
else
    echo "Compilation failed"
fi
