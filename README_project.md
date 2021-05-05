# Optimal Solutions in Rayleigh-Benard Convection: Parallel Implementation
<span style="color:red">Professor Sondak's original header:</span>

This code finds exact, steady solutions that optimize vertical heat transport in 2D Rayleigh-Benard convection. Specific
details on the approach and the numerics can be found in [[1]](#1).


## 1. Problem Description and Need for HPC
*Description of problem and the need for HPC and/or Big Data*

<span style="color:red">Pull from **Proposal Presentation** and Prof. Sondak's papers</span>
- High level overview (more detail in section #3)
- 2D Rayleigh-Benard convection
- HPC problem: compute-intensive, originally full serial


## 2. Description of Solution and Comparison to Existing Work
*Description of solution and comparison with existing work on the problem*
- Reference Prof. Sondak's paper
- I think "solution" here means 2-3 sentences about parallelized implementation
- <span style="color:red">Additional papers?</span>


## 3. Model Description in Detail
*Description of your model and/or data in detail: where did it come from, how did you acquire it, what does it mean, etc.*
- Fluid between horizontal plates with fixed temperature difference. 
- Competing buoyancy and gravitational forces.
- Convection cells, flow characterized by Rayleigh number. 
- Nusselt number as ratio of convective to conductive heat transfer.
- System modelled by Oberbeck-Boussinesq equations of motion. 
- Fourier transform along x. Chebyshev finite difference along y.
- Finer mesh at high Rayleigh, more computationally intensive
- Parallel processing aids higher Rayleigh, turbulent simulation. 


## 4. Technical Details of Parallel Design
*Technical description of the parallel application, programming models, platform and infrastructure*
### 4.1 Parallel Application
- HPC again?
- Necessary to hit levels/types from Design Presentation?
### 4.2 Programming Models
- OpenMP & MPI?
### 4.3 Platform and Infrastructure
- AWS?
- Cluster?


## 5. Source Code
*Links to repository with source code, evaluation data sets and test cases*
- Quick orientation to repo, directories, code?


## 6. Technical Description of Code
*Technical description of the software design, code baseline, dependencies, how to use the code, and system and environment needed to reproduce your tests*
### 6.1 Software Design
- Pertinent components of Design Presentation
### 6.2 Code Baseline
- Put this as 6.1? Is this talking about profiling ?--if so add the Mathcha profile image. Or skip this because it's covered in section 5?
### 6.3 Dependencies
- Anything specific to Fortran or this code base?  FFTW, OpenMP, MPI?
### 6.4 Compiling and Running the Code
- Take some/all of Prof. Sondak's description below.
### 6.5 Replicability Information
- Table? 
- Input File, AWS specs, Compiler, Libraries, (Cluster Info?)
- Need a set Ra, Nx, Ny series to run.


## 7. Performance Evaluation (aka. Results?)
*Performance evaluation (speed-up, throughput, weak and strong scaling) and discussion about overheads and optimizations done*
- Use framework from Design Presentation to address points we were concerned about, and then what actually happened.
- Have screenshots of speedups graphs here


## 8. Advanced Features
*Description of advanced features like models/platforms not explained in class, advanced functions of modules, techniques to mitigate overheads, challenging parallelization or implementation aspects...*
- FFTW?
- Fortran aspects?


## 9. Discussion and Future Work
*Final discussion about goals achieved, improvements suggested, lessons learnt, future work, interesting insights…*
### 9.1 Summary of Results --OR-- Conclusions
- list of accomplishments
### 9.2 Challenges --OR-- Lessons Learned
- Fluid flow, Fortran, FFTW...AWS > Cluster
### 9.3 Future Work
- Parallel algorithm for a solving a tridiagonal matrix

## 10. References
<a id="1">[1]</a> 
Sondak, D., Smith, L.M., Waleffe, F. (2015). 
Optimal heat transport solutions for Rayleigh-Benard convection
Journal of Fluid Mechanics, 784, 656-595.
- Add FFTW?
- Add RK_IMEX method paper?
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

