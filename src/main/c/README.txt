/**
@mainpage SLIM-curve package for exponential curve fitting of spectral lifetime data.

SLIM Curve is a curve fitting library used for Fluorescent Lifetime Imaging or FLIM and Spectral Lifetime Imaging or SLIM.  It is based on code developed by Paul Barber and the Advanced Technology Group at the Gray Institute for Radiation Oncology & Biology, University of Oxford, and used for FLIM functionality in his TRI2 (Time Resolved Imaging) software.  It is also used in the LOCI SLIM Plugin project.

For exponential lifetime fitting there are two core algorithms within SLIM Curve: The first is a triple integral method that does a very fast estimate of a single exponential lifetime component. The second is a Levenberg-Marquardt algorithm or LMA that uses an iterative, least-squares-minimization approach to generate a fit. This works with single, double and triple exponential models, as well as stretched exponential. There is also code to perform 'global' analysis over a number of signals symultaneously (e.g. over an image), where the lifetimes can be considered constant across the data set, but the amplitudes are allowed to vary for each signal. There is also a completely generic global analysis function. A third algorithm is available to perform phasor analysis.

In addition there is a non-negative linear least squares algorithm that is useful for spectral unmixing in SLIM.

The code is written in C89 compatible C and is thread safe for fitting multiple pixels concurrently.  Several files are provided as wrappers to call this library from Java code:  EcfWrapper.c and .h provide a subset of function calls used by SLIM Plugin, these may be invoked directly from Java using JNA.  In addition there is a Java CurveFitter project that provides a wrapper to the SLIM Curve code.  This invokes the C code using JNI, with loci_curvefitter_SLIMCurveFitter.c and .h.

For further details, see:
    http://loci.wisc.edu/software/slim-curve
    
If you are familiar with the program TRI2, that uses SLIM Curve, this screenshot may help you to understand the meaning of the parameters.

@image html params_in_tri2.png
@image latex params_in_tri2.jpg "How some SLIM-Curve paramters are used in TRI2." width=15cm
    
    
# Library Contents:

Directory            | Contents
---------            | --------
src	                 | source files
src/main/c           | The source files for the SLIM Curve library
src/slim-curve-cmd/c | The source files for the stand alone executable wrapper for the library
test_files           | dat and ini settings file for testing
src/main/c/doc       | API documentation (Doxygen output)



# To Build the Stand Alone Program using CMake and gcc under Linux:

Create a build folder, and cd to it

    mkdir build
    cd build

Run CMake

    cmake ../CMakeLists.txt
    
Run make

    make
	
# To Run the Stand Alone Executable
Copy the executable to the test_files folder for convenience

    cp slim-curve-cmd ../test_files
    
Run the program with the test files

    cd ../test_files
    ./slim-curve-cmd test.ini transient.dat
    
    
\author Paul Barber
\copyright Creative Commons BY-SA 2014
    
*/
