
~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████  
████  ████ ██   ██ ██    ██         ██ ██       
██ ████ ██ ███████ ██    ██     █████  ███████  
██  ██  ██ ██   ██ ██    ██         ██ ██    ██ 
██      ██ ██   ██ ██    ██    ██████   ██████         
~~~


#### GPU-based Finite difference code for DNS of Multiphase Homogenous isotropic turbulence (PFM and particles)

Main developer: A. Roccon 

Future developments:
* Conservative Allen-Cahn equation for phase-field modeling.
* Particles tracking (interpolation is done)
* FFTW for CPU debug run 

Performance (only NS)
* 64  x  64 x  64 | RTX5000 |   5 ms/timestep
* 128 x 128 x 128 | RTX5000 |  30 ms/timestep
* 256 x 256 x 256 | RTX5000 | 240 ms/timestep

#### Systems supported:
* Unix + nvfortran 

#### To run a simulation:
* go to src folder and run ./go.sh

#### Parallelization strategy
* The code is serial and exploit a single GPU, the poisson solver can be extended to use all the GPUs on the node 

#### Output and restart files.
* Files containing the Eulerian fields (u\_\*\*\*,v\_\*\*\*\*,w\_\*\*\*\*  etc.) and the particle positions (p\_\*\*\*) are stored in set_run/results

