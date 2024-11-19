gcc -o ex1 ex1.c MathTool.c
./ex1 > result.dat
diff -u result.dat answer.dat
