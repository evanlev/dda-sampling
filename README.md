# On-the-fly Adaptive k-Space Sampling Using Moment-based Spectral Analysis

Update 3/4/19: I am in the process of cleaning up the C++ version and creating cmake files. 
This was originally ported from the C version and did not fully exploit C++ constructs.

DISCLAIMER: This is research done at Stanford University. We do not 
guarantee that the code is error-free. Please contact the author if you 
find any bugs or have any suggestions or thoughts.
Evan Levine (egl@stanford.edu) 2017.

## 1. Components

There are MATLAB, C, and C++ implementations.

### matlab/
MATLAB code, test scripts, and two demos. Implementations of sampling functions
in MATLAB are provided in matlab/dda_utils.

### BART/
C code that can be used with the Berkeley Advanced Reconstruction Toolbox (BART).
If this is not set up, MATLAB code will revert to slower implementations in MATLAB. 
The following functions are provided.

##### dda_getw.c
Calculate w from sensitivity maps.

##### dda_getp.c
Calculate p from sampling pattern.

##### dda_getDeltaJ.c
Calculate Delta J from w and p.

##### dda_bc.c
Best candidate sampling and approximate best candidate samplign.

### cpp/
C++ implementation of best candidate sampling algorithms and function to compute Delta J.
The BART implementation has the same functionality and a few more algorithms. I preferred
C++ to implement the data structures used in best candidate sampling.

## 2. Installation
===============

### 2.1. Prerequisites

MATLAB. 

g++ and BART highly recommended for functions to run fast. 
For cpp version, cmake and boost are required. On macOS, run "brew install boost" with homebrew.

### 2.2. Setup

##### cpp/
t can compile on Ubuntu 14.04 with g++ 4.8.5

##### BART/

Steps to set up BART code are:

1) Install BART. Your BART directory, should have a Makefile and subdirectory src. 
   I'll refer to this directory as "bart."
2) Copy the directory dda_tools to bart/src/. 
   This is the library with sampling functions
3) Copy the C files in bart, dda_getp.c, dda_getw.c, dda_getDeltaJ.c, dda_bc.c 
   to bart/src/. 
4) Add these 4 lines to the Makefile with the other MODULES_*. It should look like
   BART/Makefile

MODULES_dda_bc = -ldda_tools 
MODULES_dda_getw = -ldda_tools 
MODULES_dda_getp = -ldda_tools 
MODULES_dda_getDeltaJ = -ldda_tools 

Modify XTARGETS in Makefile
XTARGETS += $(TBASE) $(TFLP) $(TNUM) $(TIO) $(TRECO) $(TCALIB) $(TMRI) $(TSIM) $(TSAMPLING)

5) Copy dda_tools.mk to bart/rules

6) Add the line to build_targets.mk
TSAMPLING=dda_bc dda_getp dda_getw dda_getDeltaJ

7) Compile the BART code with:
    $ make dda_getp
    $ make dda_getw
    $ make dda_getDeltaJ
    $ make dda_bc


