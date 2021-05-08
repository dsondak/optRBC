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

echo "########################"
echo "#  MPI strong scaling  #"
echo "########################"

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

rm Nu_data.txt
rm input.data
