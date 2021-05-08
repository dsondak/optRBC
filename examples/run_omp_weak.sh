cd ..
make
cp time_loop.exe examples/
cd examples/

text="7.0, 1.5585, 0.1
5.0e+03, 1, 1.1
100., 0.01
-1.0, 1.0
1600, 1200, 1
0, 0, 0
0, 0, 0"
echo -e "$text" > input.data

echo "########################"
echo "#  omp weak scaling    #"
echo "########################"

echo "Serial version running ..."
export OMP_NUM_THREADS=1
timeout --foreground 100s ./time_loop.exe
echo "2 threads running ..."
export OMP_NUM_THREADS=2
timeout --foreground 100s ./time_loop.exe
echo "3 threads running ..."
export OMP_NUM_THREADS=3
timeout --foreground 100s ./time_loop.exe
echo "4 threads running ..."
export OMP_NUM_THREADS=4
timeout --foreground 100s ./time_loop.exe
echo "5 threads running ..."
export OMP_NUM_THREADS=5
timeout --foreground 100s ./time_loop.exe
echo "6 threads running ..."
export OMP_NUM_THREADS=6
timeout --foreground 100s ./time_loop.exe
echo "7 threads running ..."
export OMP_NUM_THREADS=7
timeout --foreground 100s ./time_loop.exe
echo "8 threads running ..."
export OMP_NUM_THREADS=8
timeout --foreground 100s ./time_loop.exe
echo "done!"

rm Nu_data.txt
rm input.data
