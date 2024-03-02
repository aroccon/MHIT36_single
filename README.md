
~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████  
████  ████ ██   ██ ██    ██         ██ ██       
██ ████ ██ ███████ ██    ██     █████  ███████  
██  ██  ██ ██   ██ ██    ██         ██ ██    ██ 
██      ██ ██   ██ ██    ██    ██████   ██████         
~~~


#### GPU-based Finite difference code for DNS of Multiphase Homogenous isotropic turbulence (PFM and particles)

Main developer: A. Roccon 

Current capabiltiies:
* DNS of single-phase flow
* Tracking of tracers

Future developments:
* Conservative Allen-Cahn equation for phase-field modeling.
* Tracking of inertial particles
* Improve integration of NS (from EE to AB)
* FFTW for CPU debug run 

Performance (only NS)
* 64  x  64 x  64 | RTX5000 |   1 ms/timestep
* 128 x 128 x 128 | RTX5000 |   7 ms/timestep
* 256 x 256 x 256 | RTX5000 |  50 ms/timestep
* 384 x 384 x 384 | RTX5000 | 180 ms/timestep

#### Systems supported:
* Unix + nvfortran 

#### To run a simulation:
* go to src folder and run ./go.sh

#### Parallelization strategy
* The code is serial and exploit a single GPU, the poisson solver can be extended to use all the GPUs on the node 

#### Output files.
* Files containing the Eulerian fields (u\_\*\*\*,v\_\*\*\*\*,w\_\*\*\*\*  etc.) and the particle positions (xp\_\*\*\*) are stored in src/output

