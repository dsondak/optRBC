cd ../..
make
cp re_to_comp_test.exe examples/fft_examples
cd examples/fft_examples

echo "Starting"
./re_to_comp_test.exe > one_sided.txt
echo "Done"
python one_sided_grapher.py
