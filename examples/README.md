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

## Examples

We present the 7 examples below. 

1. OpenMP strong scaling.
2. OpenMP weak scaling. 
3. Hybrid OpenMP + MPI performance.
4. MPI strong scaling.
5. MPI weak scaling.
6. Nusselt number calculation.
7. Temperature field visualization.





