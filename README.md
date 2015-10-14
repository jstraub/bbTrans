This library provides basic algorithms and data structures to
efficiently interface with GPU data as well as several helper classes
for evaluation and logging.

### Dependencies

This code depends on the following other libraries and was tested under Ubuntu
14.04. 
- Eigen3 (3.0.5) 
- cuda 5.5 or 6.5 

The GPU kernels were tested on a Nvidia Quadro K2000M with compute
capability 3.0.

### Library
*libjsCore.so* collects all the cuda code into one shared library. The rest
of the code is in the form of header files.

