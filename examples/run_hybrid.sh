cd ..
make
cp time_loop.exe examples/
cp time_loop_MPI.exe examples/
cd examples/

text="7.0, 1.5585, 0.1
5.0e+03, 1, 1.1
0.1, 0.01
-1.0, 1.0
1280, 960, 1
0, 0, 0
0, 0, 0"
echo -e "$text" > input.data

export OMPI_MCA_btl=^openib

echo "Serial version running ..."
export OMP_NUM_THREADS=1
time ./time_loop.exe >> out
echo "2 processors 2 threads running ..."
export OMP_NUM_THREADS=2
time mpirun -np 2 ./time_loop_MPI.exe >> out
echo "4 threads running ..."
export OMP_NUM_THREADS=4
time ./time_loop.exe >> out
echo "4 processors running ..."
export OMP_NUM_THREADS=1
time mpirun -np 4 ./time_loop_MPI.exe >> out
echo "2 processors 4 threads running ..."
export OMP_NUM_THREADS=4
time mpirun -np 2 ./time_loop_MPI.exe >> out
echo "4 processors 2 threads running ..."
export OMP_NUM_THREADS=2
time mpirun -np 4 ./time_loop_MPI.exe >> out
echo "8 threads running ..."
export OMP_NUM_THREADS=8
time ./time_loop.exe >> out
echo "done!"

rm Nu_data.txt
rm input.data
