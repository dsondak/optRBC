# FFT Examples.

This is a subset of the examples that correspond solely to FFT based experiments. Similar to the entire examples set, for each example there is:

1. A shell script to execute.
2. A description in this file of what is happening.

As in the general examples README, we recommend using a t2.2xlarge AWS instance to replicate results. 

## Prerequisites
 
First clone this repository to the AWS instance with 

```
git clone https://github.com/dsondak/optRBC.git
```

Change into the optRBC directory.

```
cd optRBC
```

Checkout the `parallel_project` branch with

```
git checkout parallel_project
```

Run the `aws_setup.sh` script with 

```
./aws_setup.sh
```

## Contents.

FFT examples consists of: 

1. Multithreaded FFT benchmarking.
2. One-sided FFT benchmarking. 

## 1. Multithreaded FFT

To run and visualize using python, run: 
```
./run_multithread_fft.sh
```
This iterates over various array sizes and thread counts, with thread counts ranging from 1 to 16, and executes a discrete fourier transform. The times how long it takes to execute the dft as a function of array size and thread count when planned. 

## 2. One-sided FFT

To run and visualize using python, run: 
```
./run_fft_one_side.sh
```
This iterates over various array sizes and compares the execution time of a one-sided vs. two-sided forward and backward transformation. 