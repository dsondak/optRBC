# Examples.

This directory contains a number of examples and test cases for our parallelization of the optRBC code. The examples demonstrate performance results, 
how to visualize the temperature field, and how to calculate the nusselt number. 

Each example will be composed of 3 parts:

1. A shell script to execute.
2. A description in this file of what is happening.
3. A sample output file showing what the standard output of each shell sript should look like. 

We strongly recommend that the examples are run on a `t2.2xlarge` instance on AWS with the Ubuntu 18.04 operating system, because these are the results that are presented and documented in the
sample output files. Before we jump in, a few notes.

1. In general we have found that the OpenMP results are much easier to replicate and have much more consistent runtimes.
2. We had one strange experience where a `t2.2xlarge` instance with a public ip starting in the 100 range gave very poor performance. Even though the hardware specs looked the same. If when spinning up an instance, the public IP is in that range, consider relaunching to try to get on that starts in the 30s, e.g., `35.175.132.17`.
3. These examples were run in the US East 1 (N. Virginia) region, so we enourage the reader to use that region.
4. We have noticed that the overall performance of the OpenMP code degrades after the MPI versions are run. This is very strange and tough to debug, but for the meantime we encourage the reader to run the examples in the order that they are presented, and keep in mind that once the MPI examples are run, results might be worse than before for OpenMP. Before the MPI versions are run, the OpenMP results should be extremely stable and replicable.


## Prerequisites.

The following steps must be run in order for the examples to work. 

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

This script will require you to hit enter a few times to confirm the packages being added to the instance. Once all the packages are installed run 

```
make
```

to build the binary, and 

```
cd examples
```
to enter this directory.

## Contents.

We present the 7 examples below. 

1. OpenMP strong scaling.
2. OpenMP weak scaling. 
3. Hybrid OpenMP + MPI performance.
4. Nusselt number calculation.
5. Temperature field visualization.
6. MPI strong scaling.
7. MPI weak scaling.

The remainder of this document will step through each example to demonstrate how to reproduce the results. 

### 1. OpenMP strong scaling.

This example demonstrates strong scaling according to Amdahl's law for the OpenMP multi-threaded version of the code. The problem size is
fixed to a grid of size Nx=1280 and Ny=1080. With a timestep of 0.1, we iterate for 30 steps and time the entire execution. Since the problem
size is fixed, this is an example of strong scaling. 

To run this example, run the shell script:

```
./run_omp_strong.sh
```

The total runtime of the script with be about 5 minutes on a `t2.2xlarge` instance. The file `sample_omp_strong_results` has an example output from the 
bash script to compare to. Since the linux `time` command is used, we are interested in the `real` reported time. For example, for 8 threads,

```
8 threads running ...

real	0m14.410s
user	1m29.505s
sys	0m0.356s
```
the time we are interested in is the 14.410 seconds. The figure below shows the runtimes and scalings for various problem sizes. The problem size in the 
bash script corresponds to the black circles on the plot below. 

![omp_strong](../figs/omp_strong_scaling.png)

For a discussion on these results, please refer to the main README.

### 2. OpenMP weak scaling.

This example demonstrates weak scaling according to Gustafson's law for the OpenMP multi-threaded version of the code. In order to allow the problem size to 
grow as we increase the number of threads used, we fix the grid size to Nx=1600 and Ny=1200, but allow every version of the code to run for 100 seconds. As a
result, the serial portion of the code (initialization of the matrices) is fixed for each run, and we can measure how much work each version of the code can 
do in a fixed amount of time.

To run this example, run the shell script:

```
./run_omp_weak.sh
```

The total runtime of the script will be 800 seconds (13.3 minutes). The file `sample_omp_weak_results` has an example output from the 
bash script to compare to. There is a lot of output on this example, because we need to see how many iterations each version accomplished. Since each timestep is only 0.01, we can count the number of iterations just by looking at the final time once the next version is run. For example, for the serial version, the output looks like 

```
 Running with            1  threads
 time =    1.0000000000000000E-002 dt =    1.0000000000000000E-002
 time =    2.0000000000000000E-002 dt =    1.0000000000000000E-002
...
 time =   0.25000000000000006      dt =    1.0000000000000000E-002
 time =   0.26000000000000006      dt =    1.0000000000000000E-002
2 threads running ...
```

So in 100 seconds, it was able to do 26 iterations. Now for 2 threads, 

```
Running with            2  threads
 time =    1.0000000000000000E-002 dt =    1.0000000000000000E-002
 time =    2.0000000000000000E-002 dt =    1.0000000000000000E-002
...
 time =   0.48000000000000026      dt =    1.0000000000000000E-002
 time =   0.49000000000000027      dt =    1.0000000000000000E-002
```

49 iterations were completed in 100 seconds. This process can be continued for each of the number of threads run (1-8). This data corresponds to the 
blue squares in the figure below, which plots iteration count and weak scaling speedup for two different problem sizes.

![omp_weak](../figs/omp_weak_scaling.png)

For a discussion on these results, please refer to the main README.

### 3. Hybrid OpenMP + MPI performance.

This example evaluates the performance of the hybrid implementation. The problem size is fixed at Nx=1289 and Ny=960, and is run with a step size of 0.01 for 
10 steps. This script runs 7 different experiments:

1. Serial
2. 2 processors 2 threads
3. 1 processor  4 threads
4. 4 processors 1 thread
5. 2 processors 4 threads
6. 4 processors 2 threads
7. 1 processor  8 threads

Items 2-4 are total of 4 cores, and items 5-7 are a total of 8 cores.

To run this example, run the shell script:

```
./run_hybrid.sh
```

The total runtime of the script with be about 2 minutes on a `t2.2xlarge` instance. The file `sample_hybrid_results` has an example output from the 
bash script to compare to. Again, we are interested in the `real` time.

The figure below shows the runtimes and scalings for each of the 7 configurations described above. 

![hybrid](../figs/hybrid.png)

For a discussion on these results, please refer to the main README.

