
~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████  
████  ████ ██   ██ ██    ██         ██ ██       
██ ████ ██ ███████ ██    ██     █████  ███████  
██  ██  ██ ██   ██ ██    ██         ██ ██    ██ 
██      ██ ██   ██ ██    ██    ██████   ██████         
~~~


#### GPU-based Finite difference code for DNS of Multiphase Homogenous isotropic turbulence (PFM and particles)

Main developer: A. Roccon 

Tested on:
* Milton (1 x RTX5000)
* Leonardo (1 x A100)

Current capabiltiies:
* DNS of single-phase flow
* Tracking of tracers
* Paraview files for particles
* Phase-field metho (not validated)
* Surface tension forces (not validated)


Future developments:
* Tracking of inertial particles
* Improve integration of NS (from EE to AB)
* FFTW for CPU debug run 

Performance (NS only)
* 64  x  64 x  64 | RTX5000@milton |   1 ms/timestep
* 128 x 128 x 128 | RTX5000@milton |   7 ms/timestep
* 256 x 256 x 256 | RTX5000@milton |  50 ms/timestep
* 384 x 384 x 384 | RTX5000@milton | 180 ms/timestep
* 64  x  64 x  64 | A100@Leonardo  | 0.5 ms/timestep
* 128 x 128 x 128 | A100@Leonardo  |   1 ms/timestep
* 256 x 256 x 256 | A100@Leonardo  |   7 ms/timestep
* 384 x 384 x 384 | A100@Leonardo  |  23 ms/timestep
* 512 x 512 x 512 | A100@Leonardo  |  58 ms/timestep


#### Systems supported:
* Unix + nvfortran 

#### To run a simulation:
* go to src folder and run ./go.sh

#### Parallelization strategy
* The code is serial and exploit a single GPU, the poisson solver can be extended to use all the GPUs on the node 

#### Output files.
* Files containing the Eulerian fields (u\_\*\*\*,v\_\*\*\*\*,w\_\*\*\*\*  etc.) and the particle positions (xp\_\*\*\*) are stored in src/output

