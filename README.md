# VoxMeshSurfOpenMP
Source code for the paper: Daberdaku, S. Accelerating the computation of triangulated molecular surfaces with OpenMP. J Supercomput 75, 3426–3470 (2019). https://doi.org/10.1007/s11227-019-02803-y

## Prerequisites
* g++ >= v6.3.0 (with openmp support)
* CMake >= v3.1
* Boost C++ Libraries >= v1.55.0 (required libraries: regex system program_options thread)

## Build instructions
```bash
git clone git@github.com:sebastiandaberdaku/VoxMeshSurfOpenMP.git
cd VoxMeshSurfOpenMP
mkdir build && cd build
cmake ../src
make all
```

## Running the program
The OMP_SCHEDULE environment variable has to be set before running the program. 
Run the program with the -h (--help) parameter for more details.

Example:
OMP_SCHEDULE="guided,1" ./VoxMeshSurfOpenMP 1a8r.xyzr -r4 

