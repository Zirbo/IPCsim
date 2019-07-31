# IPCsim
Programs and scripts for Inverse Patchy Colloid simulations.

## The main simulation program...
...is found in the src directory. Inside it are a file, IPC.hpp, which contains the data structure for IPCs that is used in all the other files, and three directories:
* helpers: containing useful classes for random number generation, cell lists, pair correlation functions... that are used by all the simulation programs
* 1patch: containing the JanusIPCsim project
* 2patches: containing the IPCsim project, main tool that I am using.

## How to build
To build any of 1patch or 2patches, you need to have a c++ compiler that supports c++11, and CMake 3.5 or higher. Please follow the instructions below:
* create a new directory where the program will be built;
* navigate to the build directory that you just created, then run cmake pointing to the directory of the CMakeLists.txt of the version of the program that you want to build;
* run make in your build directory.

For example, you might do
```bash
user@machine:~/gitrepos/IPCsim$ mkdir IPCsim-build
user@machine:~/gitrepos/IPCsim$ cd IPCsim-build
user@machine:~/gitrepos/IPCsim$ cmake ../src/2patches/
user@machine:~/gitrepos/IPCsim$ make -j6
```
If when compiling you have issues with stringstream::swap, please check your g++ version, there is a bug in some versions (4.something) of the compiler: even tough c++11 is officially supported, that method is not there.

## How to run
You need the executable (IPCsim or JanusIPCsim), and in the same directory you must have
* a file called input.in, of which you find commented examples in src/1patch and src/2patches;
* a file called startingstate.xyz if you wanna resume an old simulation;
* a file called staging.in if you wanna automate the simulation, of which you find an example here.
staging.in has four columns.
The first is the temperature of the stage.
The second is a negative number if the stage is not to be repeated, and a positive number which is the maximum relative deviation of the temperature; if the average temperature in the stage differs from the specifiend one by more than that amount, the stage will be repeated.
The third is 1 if you want the simulation output (trajectory.xyt, energies, g(r) and output file) to be saved, or 0 if you want them to be deleted.
The fourth is the duration of the stage.

Call the program without parameters and the inline help will explain you how to run.

## Other directories with useful things:
* 2016version:
   The old version of the program, written in 2013. Not very efficient, don't use it. It's still here because most of the correlation functions have never been converted to the new version.
* Order-Parameters:
   A program written by either Julia Fonleitner or Gernot Pauschenwein, and refined by Günther Doppelbauer, to compute the phi-4 and phi-6 order parameters. It is in Fortran, I have no idea what it does, I just adapted it to use .xyz files and use it.
* startingstate_creators
   Python scripts that create starting configurations in .xyz formats to be used by the main program.
