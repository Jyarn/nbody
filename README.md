# NBody Simulation
NBody simulation designed to run on the Scinet Teach Cluster Supercomputer. This program splits the heavy computational load across 4 different computers and uses openmp to leverage the many cores of each computer, to speed up computation. 

To run, ensure all dependencies are installed and enter `make run` in the terminal.

## Dependencies
 - openmpi
 - openmp
