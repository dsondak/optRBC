cd ..
make
cp time_loop.exe examples/
cd examples/

text="7.0, 1.5585, 0.1
5.0e+03, 1, 1.1
5.0, 0.1
-1.0, 1.0
1600, 1275, 1
0, 0, 0
0, 0, 0"
echo -e "$text" > input.data

echo "########################"
echo "#  omp strong scaling  #"
echo "########################"

echo "Serial version running ..."
export OMP_NUM_THREADS=1
time ./time_loop.exe >> out
echo "2 threads running ..."
export OMP_NUM_THREADS=2
time ./time_loop.exe >> out
echo "3 threads running ..."
export OMP_NUM_THREADS=3
time ./time_loop.exe >> out
echo "4 threads running ..."
export OMP_NUM_THREADS=4
time ./time_loop.exe >> out
echo "5 threads running ..."
export OMP_NUM_THREADS=5
time ./time_loop.exe >> out
echo "6 threads running ..."
export OMP_NUM_THREADS=6
time ./time_loop.exe >> out
echo "7 threads running ..."
export OMP_NUM_THREADS=7
time ./time_loop.exe >> out
echo "8 threads running ..."
export OMP_NUM_THREADS=8
time ./time_loop.exe >> out
echo "done!"

rm Nu_data.txt
rm input.data
