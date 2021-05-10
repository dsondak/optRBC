cd ../..
make
cp multithreading_benchmark.exe examples/fft_examples
cd examples/fft_examples

echo "Starting"
./multithreading_benchmark.exe > multithread.txt
echo "Done"
python3 multithread_grapher.py
