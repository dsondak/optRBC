


**TODO: SECTION 5**
## 5. Technical Description of Code
*Technical description of the software design, code baseline, dependencies, how to use the code, and system and environment needed to reproduce your tests*
### 5.1 Software Design
#### 5.1.1 Code Baseline
The existing code base implemented the implicit-explicit Runge-Kutta method to solve the Oberbeck-Boussinesq partial differential equations. The code 
that we worked with was the executable generated by [time_loop.f90](./time_loop.f90). The following outline highlights the main components and subroutines 
of the program. We describe the components of the baseline code to be able to clearly describe where we made changes and added parallelization and optimization. This is mostly a reference used in the description of our OpenMP and MPI implementations, so the reader is encouraged to move ahead to the summary and profile diagram below the list if they are so inclined.

1. [time_loop.f90, line 36](./time_loop.f90#L36) - reads in the `input.data` file to determine configuration of the integration.
2. [time_loop.f90, line 237](./time_loop.f90#L237) - calls to 3 subroutines to set up the mesh: `cosine_mesh`, `dxdydz`, and `y_mesh_params`. These subroutines are defined in the [mesh_pack.f90](./mesh_pack.f90) file.
3. [time_loop.f90, line 357](./time_loop.f90#L357) - call to the `imex_rk` subroutine which is defined in the [time_integrators.f90](./time_integrators.f90) file. This is the call to the implicit-explicit RK solver. Everything up to this point is just initialization, and the main work of the program is fully contained in `imex_rk`.
4. [time_integrators.f90, line 56](./time_integrators.f90#L56) - call to the `init_bc` subroutine which is defined in the [bc_setup.f90](./bc_setup.f90) file. This subroutine initializes the boundary conditions of the matrices.
5. [time_integrators.f90, line 79](./time_integrators.f90#L79) - the beginning of the time integration `do` loop. This loop is the main work done in the program, but cannot be parallelized because each time depends on the state of the system at the previous time.
6. [time_integrators.f90, line 102](./time_integrators.f90#L102) - first call to the `calc_explict` subroutine, which is defined on [line 620](./time_integrators.f90#L620). `calc_explicit` is the main bottleneck of the code. It computes the values the `Khat` variables, which are used to update the solutions to the temperature and derivative fields. The subroutine consists of six loops: 
    * Loop 1: update Khat_phi (loop specificed by a switch statement on the stage number).
    * Loop 2: calculate the derivatives in fourier space.
    * Loop 3: executes the ffts to transform everything from fourier space to physical space (this is by far the most time consuming loop)
    * Loop 4: update the non-linear terms of the solution. 
    * Loop 5: executes ffts to bring non-linear terms back into fourier space.
    * Loop 6: update Khat_T (loop specificed by a switch statement on the stage number).
7. [time_integrators.f90, line 107](./time_integrators.f90#L107) - stage 1 `do` loop. This is the first of 3 main stage loops, each of which iterate of the columns of the variables. Each stage loop is very similar in structure, and the loops calculate the values of the `K` variables, which are used to update the solutions. 
    * [line 109](./time_integrators.f90#L109): call to `calc_vari_mod` subroutine, which is defined on [line 438](./time_integrators.f90#L438). This subroutine calculates the temperature and phi fields at the first stage by solving a tridiagonal system based on the derivatives. An example of the tridiagonal solve is at [line 518](./time_integrators.f90#L518), where `dgtsv` is called (which is a LAPACK subroutine). The `calc_vari_mod` function uses different systems depending on what stage it was called in.   
    * [line 117](./time_integrators.f90#L117) call to `calc_vi_mod` subroutine, which is defined on [line 825](./time_integrators.f90#L825). This subroutine also solves a tridiagonal system, but this time to calculate the velocity field. 
    * [line 123](./time_integrators.f90#L123) call to `update_bcs_mod` subroutine, which is defined on [line 918](./time_integrators.f90#L918). This subroutine updates the boundary conditions of the phi and uy fields.
    * [line 129](./time_integrators.f90#L129) call to `calc_implicit_mod` subroutine, which is defined on [line 602](./time_integrators.f90#L602). This subroutine calculates the currect values of Kphi and KT, which will be written to the stage specific versions of thoes variables.
    * [lines 139-141](./time_integrators.f90#L139) udpdates to the phii, Ti, and uyi fields. This updates will be used to caclulate the next explicit updates.
8. [time_integrators.f90, line 148](./time_integrators.f90#L148) - second call to `calc_explicit`. This is the exact same as the previous call, but with stage 2 specificed for Loops 1 and 6. See item number 6 above for more granular details.
9. [time_integrators.f90, line 157](./time_integrators.f90#L157) - stage 2 `do` loop. The same functions as in the stage 1 `do` loop are called, just with stage 2 parameters. See item number 7 above for more granular details.
10. [time_integrators.f90, line 192](./time_integrators.f90#L192) - third call to `calc_explicit`. This is the exact same as the previous call, but with stage 3 specificed for Loops 1 and 6. See item number 6 above for more granular details.
11. [time_integrators.f90, line 201](./time_integrators.f90#L201) - stage 3 `do` loop. The same functions as in the stage 1 `do` loop are called, just with stage 3 parameters. See item number 7 above for more granular details.
12. [time_integrators.f90, line 237](./time_integrators.f90#L237) - fourth and final call to `calc_explicit`. This is the exact same as the previous call, but with stage 4 specificed for Loops 1 and 6. See item number 6 above for more granular details.
13. [time_integrators.f90, line 245](./time_integrators.f90#L245) - the actual update of the T and phi fields, using the K and Khat variables calculated over the last 3 stages of the implicit-explicit method. 
14. [time_integrators.f90, line 256](./time_integrators.f90#L256) - `do` loop to calculate the derivative fields ux and uy at the current time step.
15. [time_integrators.f90, line 288](./time_integrators.f90#L288) - if the vtk files are being written, the `write_to_vtk` function is called. This is defined in the [write_pack.f90](./write_pack.f90) file.  

##### Summary
The above list is the nitty gritty, but overall we can describe what is happening at a high level in the `imex_rk` function. For each time step, there are four matrices that need to be updated `T, phi, ux, uy`. In order to update `T` and `phi`, the 12 variables below need to be calculated (each of which is an (Ny,Nx) matrix)

* `K1_phi, K2_phi, K3_phi`
* `K1hat_phi, K2hat_phi, K3hat_phi, K4hat_phi`
* `K1_T, K2_T, K3_T`
* `K1hat_T, K2hat_T, K3hat_T, K4hat_T`

These variables make up the pieces of the Runge-Kutta method. There are 3 stages of the code. The main stage loops (items 7,9,11 above) calculate the K variables,
and the calls to `calc_explicit` subroutine (items 6,8,10,12 above) calculate the Khat variables.

#### Profile
Now that we have a description of what the baseline code actually does, we can show a profile we created to highlight the performance of each component of the 
code. 

![profile](./figs/profile.png)

In the figure above we breakdown the components of the `imex_rk` function. Each rectangle represents a chunk of code or a subroutine, and the parenthetical numbers indicate the proportion of the time at that depth that was spent on that piece. For example, at a depth of 1, we break the runtime down into four pieces: (1) stage 1, (2) stage 2, (3) stage 3, (4) update. These portions take up 41%, 26%, 27%, and 6% of the runtime respecitvely, which totals to 100%. From this figure it is clear that the four `calc_explicit` calls consume a majority of the total runtime, which is why we break down the `calc_explicit` subroutine into 6 loops to highlight that loops 3 and 5 consume a majority of that runtime. These two loops execute the fast fourier transforms. 

### 5.3 Dependencies
- Anything specific to Fortran or this code base?  FFTW, OpenMP, MPI?
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

Future work for this project would include switching to a parallel tridiagonal solver algorithm, as well as fully debugging the Nu inconsistencies between one-sided and two-sided FFT. Moving to a parallel version of the tridiagonal solve would mitigate a bottleneck in the MPI code version, in which all of the worker nodes send their (**TODO**: which info?) to the master node. The master node then executes the tridiagonal solve and sends messages to all of the worker nodes. Having to send so many messages to the master node introduces more overhead, and having all worker nodes wait for the master node to receive all the messages would likely introduce additional idle time. Both idle time and overhead would be mitigated by using a parallel version of this tridiagonal solve.  

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

