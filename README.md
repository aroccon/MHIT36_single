
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
* DNS of single-phase flow (validated in turbulence)
* Tracking of tracers
* Paraview files for particles
* Phase-field method (not yet validated)
* Surface tension forces (not yet validated)

Future developments:
* Tracking of inertial particles
* Improve integration of NS (from EE to AB2)
* FFTW for CPU debug run 

Performance (NS only)
* 128 x 128 x 128 | RTX5000@milton |  14 ms/timestep
* 128 x 128 x 128 | A100@Leonardo  |   5 ms/timestep
* 256 x 256 x 256 | A100@Leonardo  |  40 ms/timestep

Max resolution tested:
* 256 x 256 x 256 | RTX5000@milton - 16 GB VRAM
* 768 x 768 x 768 | A100@Leonardo - 64 GB VRAm

#### Systems supported:
* Unix + nvfortran 

#### To run a simulation:
* go to src folder and run ./go.sh

#### Parallelization strategy
* The code is serial and exploit a single GPU (GPU-resident)

#### Output files.
* Files containing the Eulerian fields (u\_\*\*\*,v\_\*\*\*\*,w\_\*\*\*\*  etc.) and the particle positions (xp\_\*\*\*) are stored in src/output

