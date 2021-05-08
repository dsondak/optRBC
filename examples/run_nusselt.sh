cd ..
make
cp time_loop.exe examples/
cd examples/

echo "########################"
echo "#    nusselt number    #"
echo "########################"

export OMP_NUM_THREADS=8

text="7.0, 1.5585, 0.1
5.0e+04, 1, 1.1
100.0, 0.1
-1.0, 1.0
64, 48, 1
0, 0, 0
0, 0, 1"
echo "$text" > input.data

./time_loop.exe

mv Nu_data.txt nu1.txt

text="7.0, 1.5585, 0.1
5.0e+03, 1, 1.1
100.0, 0.1
-1.0, 1.0
64, 48, 1
0, 0, 0
0, 0, 1"
echo "$text" > input.data

./time_loop.exe

mv Nu_data.txt nu2.txt

text="7.0, 1.5585, 0.1
5.0e+02, 1, 1.1
100.0, 0.1
-1.0, 1.0
64, 48, 1
0, 0, 0
0, 0, 1"
echo "$text" > input.data

./time_loop.exe

mv Nu_data.txt nu3.txt

rm input.data
