$nodes = 8;
$matrixsize = 8;
mpiexec -n $nodes ConsoleApplication1.exe "matrix.dat" "matrix2.dat" $matrixsize 