echo "Serial version running ..."
export OMP_NUM_THREADS=1
timeout 100s ./time_loop.exe
echo "2 threads running ..."
export OMP_NUM_THREADS=2
timeout 100s ./time_loop.exe
echo "3 threads running ..."
export OMP_NUM_THREADS=3
timeout 100s ./time_loop.exe
echo "4 threads running ..."
export OMP_NUM_THREADS=4
timeout 100s ./time_loop.exe
echo "5 threads running ..."
export OMP_NUM_THREADS=5
timeout 100s ./time_loop.exe
echo "6 threads running ..."
export OMP_NUM_THREADS=6
timeout 100s ./time_loop.exe
echo "7 threads running ..."
export OMP_NUM_THREADS=7
timeout 100s ./time_loop.exe
echo "8 threads running ..."
export OMP_NUM_THREADS=8
timeout 100s ./time_loop.exe
echo "done!"
