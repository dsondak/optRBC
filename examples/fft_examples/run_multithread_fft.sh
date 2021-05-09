sudo apt install python3-pip
pip3 install numpy matplotlib pandas

cd ../..
make
cp multithreading_benchmark.exe examples/fft_examples
cd examples/fft_examples

echo "Starting"
./multithreading_benchmark.exe > multithread.txt
echo "Done"
python3 multithread_grapher.py
