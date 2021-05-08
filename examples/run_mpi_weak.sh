cd ..
make
cp time_loop.exe examples/
cp time_loop_MPI.exe examples/
cd examples/

text="7.0, 1.5585, 0.1
5.0e+03, 1, 1.1
100., 0.01
-1.0, 1.0
1280, 960, 1
0, 0, 0
0, 0, 0"
echo -e "$text" > input.data

export OMPI_MCA_btl=^openib

echo "########################"
echo "#  mpi weak scaling    #"
echo "########################"

echo "Serial version running ..."
export OMP_NUM_THREADS=1
timeout --foreground 30s ./time_loop.exe
echo "2  processes running ..."
mpirun --timeout 30 -np 2 ./time_loop_MPI.exe
echo "3  processes running ..."
mpirun --timeout 30 -np 3 ./time_loop_MPI.exe
echo "4  processes running ..."
mpirun --timeout 30 -np 4 ./time_loop_MPI.exe
echo "5  processes running ..."
mpirun --timeout 30 -np 5 ./time_loop_MPI.exe
echo "6  processes running ..."
mpirun --timeout 30 -np 6 ./time_loop_MPI.exe
echo "done!"

rm Nu_data.txt
rm input.data
