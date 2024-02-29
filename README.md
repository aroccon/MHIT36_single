
~~~text
          _____                    _____                    _____                _____          
         /\    \                  /\    \                  /\    \              /\    \         
        /::\____\                /::\____\                /::\    \            /::\    \        
       /::::|   |               /:::/    /                \:::\    \           \:::\    \       
      /:::::|   |              /:::/    /                  \:::\    \           \:::\    \      
     /::::::|   |             /:::/    /                    \:::\    \           \:::\    \     
    /:::/|::|   |            /:::/____/                      \:::\    \           \:::\    \    
   /:::/ |::|   |           /::::\    \                      /::::\    \          /::::\    \   
  /:::/  |::|___|______    /::::::\    \   _____    ____    /::::::\    \        /::::::\    \  
 /:::/   |::::::::\    \  /:::/\:::\    \ /\    \  /\   \  /:::/\:::\    \      /:::/\:::\    \ 
/:::/    |:::::::::\____\/:::/  \:::\    /::\____\/::\   \/:::/  \:::\____\    /:::/  \:::\____\
\::/    / ~~~~~/:::/    /\::/    \:::\  /:::/    /\:::\  /:::/    \::/    /   /:::/    \::/    /
 \/____/      /:::/    /  \/____/ \:::\/:::/    /  \:::\/:::/    / \/____/   /:::/    / \/____/ 
             /:::/    /            \::::::/    /    \::::::/    /           /:::/    /          
            /:::/    /              \::::/    /      \::::/____/           /:::/    /           
           /:::/    /               /:::/    /        \:::\    \           \::/    /            
          /:::/    /               /:::/    /          \:::\    \           \/____/             
         /:::/    /               /:::/    /            \:::\    \                              
        /:::/    /               /:::/    /              \:::\____\                             
        \::/    /                \::/    /                \::/    /                             
         \/____/                  \/____/                  \/____/                              
                                                                                 
~~~


#### GPU-based Finite difference code for DNS of Multiphase Homogenous isotropic turbulence (PFM and particles)

Main developer: A. Roccon 

Future developments:
* Allen-Cahn equation for phase-field modeling.
* Particles tracking

#### Systems supported:
* UNIX (ubuntu 18.04) + nvfortran 

#### To run a simulation:
* Execute ./compile.sh

#### Parallelization strategy
* The code is serial, it exploit a single GPU, the poisson solver can be extended to use all the GPUs on the node 

#### Output and restart files.
* Files containing the Eulerian fields (u\_\*\*\*,v\_\*\*\*\*,w\_\*\*\*\*,sc\_\*\*\*,  etc.) and the particle positions (p\_\*\*\*) are stored in set_run/results


#### Visualization with Paraview
Two possible cases: single-phase and particles-laden flow.
* Single-phase flow:
  - Go to set_run/results/paraview_fields and run go.sh, paraview files are generated in the output folder.
* Particles-laden flow:
  - Go to set_run/results/paraview_fields and run go.sh, paraview files are generated in the output folder.
  - Go to set_run/results/paraview_particles and run go.sh, paraview files are generated in the output folder.


#### Validation database:
https://torroja.dmt.upm.es/turbdata/agard/chapter3/HOM03/


