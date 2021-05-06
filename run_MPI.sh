time ./time_loop.exe > out
time mpirun -np 2 ./time_loop_MPI.exe > out
time mpirun -np 3 ./time_loop_MPI.exe > out
time mpirun -np 4 ./time_loop_MPI.exe > out
time mpirun -np 5 ./time_loop_MPI.exe > out
time mpirun -np 6 ./time_loop_MPI.exe > out
