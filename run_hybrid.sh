echo "Serial version running (expect 25 seconds) ..."
export OMP_NUM_THREADS=1
time ./time_loop.exe > out
echo "2 processors 2 threads running ..."
export OMP_NUM_THREADS=2
time mpirun -np 2 ./time_loop_MPI.exe > out
echo "4 threads running ..."
export OMP_NUM_THREADS=4
time ./time_loop.exe > out
echo "4 processors running ..."
export OMP_NUM_THREADS=1
time mpirun -np 4 ./time_loop_MPI.exe > out
echo "2 processors 4 threads running ..."
export OMP_NUM_THREADS=4
time mpirun -np 2 ./time_loop_MPI.exe > out
echo "4 processors 2 threads running ..."
export OMP_NUM_THREADS=2
time mpirun -np 4 ./time_loop_MPI.exe > out
echo "8 threads running ..."
export OMP_NUM_THREADS=8
time ./time_loop.exe > out
echo "done!"
