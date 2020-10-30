# Optimal Solutions in Rayleigh-Benard Convection

This code solves the 2D Boussinesq equations between two infinite plates. The wall-normal direction is discretized using
non-uniform finite differences. The x-direction uses a Fourier basis.

## Code Overview

The flowmap algorithm is used to find steady solutions to the 2D Boussinesq equations.  In addition to this, a simple
optimization is performed to find the maximum Nu at a given Ra by varying the "size of the box"; that is, by restricting the
size of permitted structures.

The optimization is done by simply looping over various values of $\alpha$, sensing when a maximum has been crossed in Nu,
and fitting a parabola around that maximum.

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
* `vtkdata/`: Solution fields are written out in VTK format at each Ra for solutions optimizing heat transport
* `uy`, `temperature`: Vertical velocity and temperature fields in binary format. Used for restarts. Note that `ux` can be
computed from `uy` in 2D from the continuity equation.
