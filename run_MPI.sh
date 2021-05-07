export OMP_NUM_THREADS=1
echo "Serial version running ..."
time ./time_loop.exe > out
echo "2 processors running ..."
time mpirun -np 2 ./time_loop_MPI.exe > out
echo "3 processors running ..."
time mpirun -np 3 ./time_loop_MPI.exe > out
echo "4 processors running ..."
time mpirun -np 4 ./time_loop_MPI.exe > out
echo "5 processors running ..."
time mpirun -np 5 ./time_loop_MPI.exe > out
echo "6 processors running ..."
time mpirun -np 6 ./time_loop_MPI.exe > out
echo "done!"
