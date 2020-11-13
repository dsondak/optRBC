# Optimal Solutions in Rayleigh-Benard Convection

This code finds exact, steady solutions that optimize vertical heat transport in 2D Rayleigh-Benard convection. Specific
details on the approach and the numerics can be found in [[1]](#1).

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
* `vtkdata/`: Solution fields are written out in VTK format at each Ra for solutions optimizing heat transport
* `uy`, `temperature`: Vertical velocity and temperature fields in binary format. Used for restarts. Note that `ux` can be
computed from `uy` in 2D from the continuity equation.

## Contributing to the Project
Contributions to the code are very welcome. To contribute, please do the following:
1. Create a new issue
2. Move that issue to "In Progress" on the Kanban board
3. Create a new branch
    1. To create a branch locally do `git checkout -b <new_branch_name>`
    2. To push this branch to the main repo do `git push --set-upstream origin <new_branch_name>`
        * This is only necessary the first time you create a branch.
4. Make changes to the code on your branch and push your changes when they're ready 
    1. `git commit` locally
    2. `git push` to push your new changes
    3. **Recommendation:** Only push changes that compile
5. When you're ready, make PR to `main`
    * Assign a reviewer
6. After making the PR, link the issue with that PR
7. Continue making changes on that branch until they are approved to be merged to `main`

## References
<a id="1">[1]</a> 
Sondak, D., Smith, L.M., Waleffe, F. (2015). 
Optimal heat transport solutions for Rayleigh-Benard convection
Journal of Fluid Mechanics, 784, 656-595.
