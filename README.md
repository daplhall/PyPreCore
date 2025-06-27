# A post-processing dust simulator for gas models    
PyPreCore was written as part of my thesis which can be downloaded [here](https://daplhall.github.io/pdf/DanielHallThesis.pdf). Its a post-processing dust simulation to
run with a gas background model which can be either Lagrangian and Eulerian.

# Usage
The main class of the program is `Particles` which represents the dust particles.
One of its inputs is an environment, which is either a `DataSphere` or `BE_sphere` (Bonnor-ebert sphere). The particles and its environment interact through the `drag` method.  
  