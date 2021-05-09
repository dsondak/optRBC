sudo apt install python3-pip
pip3 install numpy matplotlib pandas meshio

cd ..
make
cp time_loop.exe examples/
cp time_loop_MPI.exe examples/
cd examples/

export OMPI_MCA_btl=^openib

mkdir vtkdata

text="7.0, 1.5585, 0.1
5.0e+04, 1, 1.1
5.0, 0.1
-1.0, 1.0
64, 48, 1
0, 0, 0
1, 0, 0"
echo "$text" > input.data

echo "########################"
echo "#      temp viz        #"
echo "########################"

./time_loop.exe

mpirun -np 4 ./time_loop_MPI.exe

rm Nu_data.txt
rm input.data
