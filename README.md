SLIM Curve is a curve fitting library used for Fluorescent Lifetime Imaging or
FLIM and Spectral Lifetime Imaging or SLIM.  It is developed by Paul Barber and
the Advanced Technology Group at the [Cancer Research UK and Medical Research
Council Oxford Institute for Radiation Oncology] (http://www.rob.ox.ac.uk/),
University of Oxford, as well as the [Laboratory for Optical and Computational
Instrumentation](http://loci.wisc.edu/) at the University of Wisconsin-Madison.
SLIM Curve is used for FLIM functionality in the Advanced Technology
Group's [Time Resolved Imaging](https://www.assembla.com/spaces/ATD_TRI/wiki)
(TRI2) software, as well as in the [SLIM Curve plugin for
ImageJ](http://fiji.sc/SLIM_Curve).

For exponential lifetime fitting there are two core algorithms within SLIM
Curve:

1. A triple integral method that does a very fast estimate of a single
   exponential lifetime component.
2. The second is a Levenberg-Marquardt algorithm or LMA that uses an iterative,
   least-squares-minimization approach to generate a fit. This works with
   single, double and triple exponential models, as well as stretched
   exponential.

There is also code to perform 'global' analysis over a number of signals
simultaneously (e.g. over an image), where the lifetimes can be considered
constant across the data set, but the amplitudes are allowed to vary for each
signal. There is also a completely generic global analysis function. A third
algorithm is available to perform phasor analysis.

In addition there is a non-negative linear least squares algorithm that is
useful for spectral unmixing in SLIM.

The code is written in C89 compatible C and is thread safe for fitting multiple
pixels concurrently. Several files are provided as wrappers to call this
library from Java code: `EcfWrapper.c` and `.h` provide a subset of function
calls used by the ImageJ plugin; these may be invoked directly from Java using
JNA. In addition there is a Java CurveFitter project that provides a wrapper to
the SLIM Curve code. This invokes the C code using JNI, with
`loci_curvefitter_SLIMCurveFitter.c` and `.h`.

## See also

* [SLIM Curve wiki](https://github.com/slim-curve/slim-curve/wiki)
* [SLIM Curve web site](https://slim-curve.github.io/)
* [SLIM Curve Doxygen docs](http://code.imagej.net/slim-curve/html/)

## Directory contents

* `src` - source files
* `src/main/c` - The source files for the SLIM Curve library
* `src/slim-curve-cmd/c` - The source files for the standalone executable
  wrapper for the library
* `src/slim-curve-cmd/c/doc` - API documentaion
* `test_files` - dat and ini settings file for testing

## Building the source

To build the standalone program using CMake and gcc under Linux:

1.  Create a build folder, and cd to it

    ```
    mkdir build
    cd build
    ```

2.  Run CMake

    ```
    cmake ../CMakeLists.txt
    ```

3.  Run make

    ```
    make
    ```

## Running the standalone executable

1.  Copy the executable to the `test_files` folder for convenience

    ```
    cp slim-curve-cmd ../test_files
    ```

2.  Run the program with the test files

    ```
    cd ../test_files
    ./slim-curve-cmd test.ini transient.dat
    ```
