cd ../..
make
cp re_to_comp_test.exe examples/fft_examples
cd examples/fft_examples

./re_to_comp_test.exe > one_sided.txt
python one_sided_grapher.py
