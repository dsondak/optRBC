# Optimal Solutions in Rayleigh-Benard Convection: Parallel Implementation

This `parallel_project` branch captures the results of a semester long project in Harvard's CS205: Computing Foundations for Computational Science. The primary goal was to modify the time integration step of the code to compute calculations in parallel.  The code base and inspiration for the project came from Sondak et al., of the paper below [[1]](#1). Team Members: Katrina Gonzalez, Michael Neuder, Jack Scudder


## 1. Problem Description and Need for HPC
*Description of problem and the need for HPC and/or Big Data*

### 1.1 Model Description in Detail
*Description of your model and/or data in detail: where did it come from, how did you acquire it, what does it mean, etc.*

Rayleigh-Benard convection (RBC) is phenomena that takes place when a liquid is placed between two (approximately) infinite plates held at some temperature difference. When the bottom plate is kept at a higher temperature than the top plate, this creates competing buoyancy and gravitational effects.  A simple and complex example are below.

![rbc_image](figs/rbc_simple.png)

The temperature flow of this liquid can be described by the Oberbeck-Boussinesq partial differential equations. These PDEs can be solved at each location and time step to numerically compute the temperature field and its evolution. The non-dimensional parameters we are primarily focused on are listed with descriptions below:
- Rayleigh Number (Ra) is a measure of convection. 
- Prandtl Number (Pr) is a measure of viscosity (eg. Pr=7 for water, Pr=98.31 for engine oil). 
- Nusselt Number (Nu) is a ratio of convective to conductive heat transfer. 

![nondim_params](figs/ra_pr_nu_eqns.png)
 \* adapted from a presentation based on the paper in [[1]](#1).

The low Rayleigh number regime in RBC has been studied extensively; therefore, this code was written to simulate problems with high Rayleigh numbers.  Though Rayleigh-Benard convection is a physical and therefore 3D problem, this particular code implements 2D Rayleigh-Benard convection.  This simplification reduces the complexity of the algorithm, allowing for exploration of the ratio between Nu and Ra at high values of Ra up to 10^9. 

### 1.2 Need for HPC
At high Ra numbers the discretized mesh is finer, and therefore more computationally intensive to run the simulation.  In order to simulate high Ra numbers more quickly, we sought to parallelize the existing code using OpenMP and MPI, and we explored Fast Fourier Transform (FFT) methods to improve algorithmic efficiency.  (*Note*: our particular FFT implementation is FFTW3 [[2]](#2).  We use FFT and FFTW interchangeably to reference the same calculations.)  The serial version of the current algorithm has a time complexity of O(N^3 * log(N)), leading to exceedingly long computation times at large problem sizes.  Since this is a compute-intensive code designed to solve complex PDEs, it is a High Performance Computing (HPC) problem.


## 2. Description of Solution and Comparison to Existing Work
*Description of solution and comparison with existing work on the problem*

Using the Fortran code base from Sondak et al. [[1]](#1), we initially profiled the serial code to determine the primary bottlenecks of the `time_integrators.f90` portion.  In each time step, the code performs updates using an implicit-explicit Runge-Kutta method detailed in [[3]](#3).  In the code, the `imex_rk()` subroutine computes 8 Fast Fourier Transforms per time step, each of which costs O(N * log(N)) computational time.  Details of the code profiling are shown below.

![serial_profile](figs/serial_profile.png)

Values in parenthesis are the percentage of that parent subroutine's runtime spent on the boxed subroutine.  The subroutine `calc_explicit()` performs the FFT calculations and is the clear bottleneck of the code.  Looking at this profile helped us identify our three-step approach to parallelization: (1) use OpenMP to parallelize nested loops at each time step, (2) use MPI to distribute calculations across multiple nodes, and (3) explore FFTW modifications to leverage multithreading capabilities or reduce computation time by using a one-sided FFT.

As previously stated, our primary project comparison is to the base code from Sondak et al. [[1]](#1).  The endstate of our code is the same as the original: to calculate optimal solutions in RBC.  We merely sought to make the process more efficient.  Without doing a full literary review, one particular paper by Clarke et al. (2020) [[4]](#4) sought to parallelize the same problem set.  They applied the same discretization using the Boussinesq approximation to the Navier-Stokes equations for a 2D RBC simulation, but used the Parareal algorithm.  They claimed to achieve speedups of up to 2.4 "with around 20 processors \[...\] for Rayleigh numbers as high as 10^6."  Their code wasn't shared; however, so we were not able to directly compare their implementation with ours.  Additionally, it was both helpful and educational to experiment with Fortran without having to create this robust simulation from scratch.


## 3. Technical Details of Parallel Design
*Technical description of the parallel application, programming models, platform and infrastructure*

### 3.1 Parallel Application
As stated in Section 1, this is a High Performance Computing problem.  The original implementation was completely serial, and with a time complexity of O(N^3 * log(N)), it requires a large amount of computing power to achieve its goal of simulating RBC at high Ra.  Our design uses both fine-grained (OpenMP) parallelism at the loop and procedure level, and coarse-grained parallelism (MPI) at the task level discretization at scale.  The types of parallelism match those levels, with OpenMP allowing us to employ function/control parallelism, and MPI allowing us to apply data parallelism.  In our MPI implementation, we make use of the Single Program-Multiple Data (SPMD) as we pass large amount of data between nodes.

### 3.2 Programming Models
Our design uses both shared memory parallelism (OpenMP) and distributed memory parallelism (MPI).  While we originally intended to create a fully hybrid model, we realized that an important part of the serial code relies on the Basic Linear Algebra Subprograms (BLAS) routine that solves tridiagonal matrices. Further information is provided in Section 5.1, but our final model allows users to run an OpenMP version of the code and an MPI version of the code separately to retain some speedups that were specifically seen by using OpenMP.

Along with introducing different levels of parallelism, we also sought to introduce some speedup with enhancements like multithreaded FFT and one-sided FFT. To explore multithreaded FFT (which is included in the FFTW library), we wrote a test script in which we compared the timing of multithreaded FFT against single threaded FFT. We found that over both thread and array size, there was no meaningful speedup observed.

With respect to one-sided FFT, we experimented with more strongly typing certain arrays as real (as physical fields are real) and using Hermitian conjugate symmetry (when computing the Fourier Transform of a size N array, only N/2 + 1 of those elements are not redundant). In particular, this change involved modifying array types, array sizes, and certain loops meant to be executed in Fourier space. Subroutines swapped out include planning and execution, with the strictly real to complex or complex to real transformations included in the FFTW library. Due to this symmetry, we expected a theoretical speedup of 2x for each Fourier transform (originally N*log(N) ), along with additional 2x speedups for Nx loops in Fourier space. 

### 3.3 Platform and Infrastructure
We used an Amazon Web Services (AWS) t2.2xlarge instance with the Ubuntu 18.04 operating system for our evaluation and performance criteria.  Further details are included in the [examples directory](https://github.com/dsondak/optRBC/tree/parallel_project/examples) README.

We planned to use to the FAS RC academic cluster to run our experiments, but ran into several challenges, two of which are captured here.  First, the results were inconsistent depending on the day of the week and the time of day.  Despite using the Slurm commands introduced on the FAS RC documentation and in class, the results were still inconsistent.  This was likely due to a high degree of device sharing between jobs, introducing variability with respect to how much of the machine we could use.  Second, the performance of our code on AWS was consistent and also yielded significantly higher speedups.  Running the OpenMP version of our code, at Nx= 1280, Ny=1080, Ra=5000 yielded a 6x speedup on an AWS instance with 8 threads.  Running the same version on the cluster yielded a ~4.5x speedup.  These results caused us to switch focus to AWS and abandon using the cluster for performance results.


## 4. Source Code
*Links to repository with source code, evaluation data sets and test cases*

<span style="color:red">**TODO: MIKE**: Please confirm accuracy</span>

Included on this `parallel_project` branch are the OpenMP and MPI implementations of the time integration step listed below:
- OpenMP:  `time_integrators.f90` and driver code `time_loop.f90`
- MPI: `time_integrators_MPI.f90` and driver code `time_loop.f90`

Modifications were made to several other files, but the majority of changes are found in those files.  FFTW modifications for one-sided computations are found in https://github.com/dsondak/optRBC/tree/real_to_comp. This particular branch and set of code changes are quarantined since one-sided FFT is not fully functional and required changes to multiple sets of files. A benchmarking code for theoretical speedup of one-sided FFT over different array size can be found <span style="color:red">**TODO: KATRINA**: finish this </span>. A benchmarking code for speedup of multithreaded FFT over several array sizes and thread counts can be found <span style="color:red">**TODO: KATRINA**: finish this </span>. 

Detailed examples and performance evaluation test cases can be found in the [examples directory](https://github.com/dsondak/optRBC/tree/parallel_project/examples).  The `README` in that directory gives stepwise instructions to compile, run, and replicate our results. 


**TODO: SECTION 5**
## 5. Technical Description of Code
*Technical description of the software design, code baseline, dependencies, how to use the code, and system and environment needed to reproduce your tests*
### 5.1 Software Design
- Pertinent components of Design Presentation
### 5.2 Code Baseline
- Put this as 6.1? Is this talking about profiling ?--if so add the Mathcha profile image. Or skip this because it's covered in section 2?
### 5.3 Dependencies
- Anything specific to Fortran or this code base?  FFTW, OpenMP, MPI? -> I don't think so

The primary dependencies of this code are on [FFTW](http://www.fftw.org/), OpenMP, and MPI. Paraview is also recommended to visualize temperature and velocity field evolution, but is not required. 

### 5.4 Compiling and Running the Code
- Take some/all of Prof. Sondak's description below.
### 5.5 Replicability Information
- Table? 
- Input File, AWS specs, Compiler, Libraries, (Cluster Info?)
- Need a set Ra, Nx, Ny series to run.

**TODO: SECTION 6**
## 6. Performance Evaluation (aka. Results?)
*Performance evaluation (speed-up, throughput, weak and strong scaling) and discussion about overheads and optimizations done*
- Use framework from Design Presentation to address points we were concerned about, and then what actually happened.
- Have screenshots of speedups graphs here with brief explanations/insights


## 7. Advanced Features
*Description of advanced features like models/platforms not explained in class, advanced functions of modules, techniques to mitigate overheads, challenging parallelization or implementation aspects...*

### 7.1 Programming in Fortran

Although Fortran is used widely used in the computational sciences for high performance computing, it was also a language not covered in the course. There were also subtleties when using Fortran with OpenMP, which included injecting variables into subroutine calls even though they were globally scoped. 

### 7.2 Working in a mathematically rich, existing codebase

We worked in an existing, extensive codebase that relied on relatively complex mathematical methods like FFT and Runge–Kutta methods to solve PDEs. This required us to study the existing code as well as its mathematical methods. One example of the need to not only understand the code but also the math implemented was revealed when we attempted to move to a one-sided FFT. Because this method relied on the symmetry of Fourier transforms of real data, it required allocating less space for complex arrays vs. real arrays. Since we also have multiple loops in Fourier space, this approach also required changing those bounds since we remove the redundant data that standard FFTs would have kept. 

**TODO: SECTION 8**
## 8. Discussion and Future Work
*Final discussion about goals achieved, improvements suggested, lessons learned, future work, interesting insights…*
### 8.1 Summary of Results --OR-- Conclusions
- list of accomplishments
### 8.2 Challenges --OR-- Lessons Learned
- Fluid flow, Fortran, FFTW...AWS > Cluster

Some challenges that we encountered with this particular project include: 
- Fortran complexities in conjunction with OpenMP, where we needed to inject variables into subroutine calls rather than relying on the global scope of those variables. 
- Small inconsistencies in values (such as Nu) between runs and whether or not those inconsistencies are due to machine precision or actual code incorrectness. This was particularly a problem when debugging the one-sided FFT approach, as it was difficult to check during intermediate steps whether or not Nu was significantly incorrect.
- Inconsistencies between runs using the Academic Cluster vs. AWS. Though it would have been appealing to run the entire project using the Academic Cluster from both a cost and (tightly coupled computers?) perspective, we ultimately chose to benchmark using AWS because it offered a more significant, more consistent level of speedup. 

### 8.3 Future Work

Future work for this project would include switching to a parallel tridiagonal solver algorithm, as well as fully debugging the Nu inconsistencies between one-sided and two-sided FFT. Moving to a parallel version of the tridiagonal solve would mitigate a bottleneck in the MPI code version, in which all of the worker nodes send their (**TODO**: which info?) to the master node. The master node then executes the tridiagonal solve and sends messages to all of the worker nodes. Having to send so many messages to the master node introduces more overhead, and having all worker nodes wait for the master node to receive all the messages introduces additional idle time. Both idle time and overhead would be mitigated by using a parallel version of this tridiagonal solve.  

Though we experimented with implementing one-sided FFT, we were ultimately unable to get it fully working. However, it would be beneficial to get this implementation working since it is compatible with both of the parallel versions (as it is a form of algorithmic speedup) and those speedups would stack multiplicatively. Methods to debug this implementation would include calculating Nu of various arrays as a kind of checksum to ensure that results are consistent between one-sided and two-sided FFT at each step in computation. 

## 9. References
<a id="1">[1]</a> 
Sondak, D., Smith, L.M., Waleffe, F. (2015). 
"Optimal heat transport solutions for Rayleigh-Benard convection."
Journal of Fluid Mechanics, 784, 656-595.

<a id="2">[2]</a> 
M. Frigo and S. G. Johnson. (2005).
"The Design and Implementation of FFTW3." 
Proceedings of the IEEE, vol. 93, no. 2, 216-231.

<a id="3">[3]</a> 
Uri M. Ascher, Steven J. Ruuth, and Brian Wetton. (1993). 
"Implicit-Explicit Methods for Time-Dependent PDE's." 
Technical Report. University of British Columbia, CAN.

<a id="4">[4]</a>
Clarke, A., Davies, C., Ruprecht, D. et al. (2020)
"Performance of parallel-in-time integration for Rayleigh Bénard convection."
Computing and Visualization in Science, 23, 10.

**TODO**
- Parallel tridiagonal solver resources shared by Prof. Sondak?
- Others?


<span style="color:red">Retained from original README:</span>

## Code Overview

The flowmap algorithm is used to find steady solutions to the 2D Boussinesq equations. The Nusselt number is computed at each
steady solution for a given $Ra$ and $Pr$. At each $(Ra, Pr)$ pair, steady solutions are computed in different box sizes,
$L$. The size of the domain (box) in $x$ is given by $L=2\pi/\alpha$ where $\alpha$ is the wavenumber that sets the scale of
the solution in $x$. The solution in box size $L$ that maximizes the Nusselt number is referred to as the optimal solution. 

### Input File Structure
| Line #   | Column 1      | Column 2                  | Column 3                          |
| :------: | :--------:    | :--------:                | :--------:                        |
| 1        | $Pr$          | Initial $\alpha$          | $\alpha$ step                     |
| 2        | Initial $Ra$  | # of $Ra$s to do          | $Ra$ increment multiplier         |
| 3        | Flow map time | Time step                 | Leave blank                       |
| 4        | $y$ at bottom | $y$ at top                | Leave blank                       |
| 5        | Nx            | Ny                        | Nz                                |
| 6        | x-refinement  | y-refinement              | z-refinement                      |
| 7        | Save to VTK   | Save binary restart files | Calculate actual optimal solution |

### Running the code
Run `make` to generate an executable. The basic executable will be called `Ra_loop.exe`. Run the executable from the command
line by typing `./Ra_loop.exe`. The code will automatically search for the input file `input.data` and read in any
parameters.

Depending on how the input file is set up, the code will loop over different values of Ra. For each Ra, it will perform an
optimization over the box size to find which box size maximizes the heat transport at the Ra. This optimization is done by
fitting a parabola to three points. The code continues until a maximum Nu is found.

The results are written out in various formats:
* `Nu.txt`: Contains the Ra, Nu, optimal wavenumber, and optimal box size
* `vtkdata/`: Solution fields are written out in VTK format at each Ra for solutions optimizing heat transport (Note: Need to create vtkdata directory before running code)
* `uy`, `temperature`: Vertical velocity and temperature fields in binary format. Used for restarts. Note that `ux` can be
computed from `uy` in 2D from the continuity equation.

