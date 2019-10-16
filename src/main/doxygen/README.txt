/**
@mainpage FLIMLib package for exponential curve fitting of fluorescence lifetime data.

FLIMLib is a curve fitting library used for Fluorescent Lifetime Imaging or FLIM.  It is based on code developed by Paul Barber and the Advanced Technology Group at the Cancer Research UK and Medical Research Council Oxford Institute for Radiation Oncology, University of Oxford, and used for FLIM functionality in their TRI2 (Time Resolved Imaging) software.  It is also used in the FLIMJ plugin for ImageJ.

For exponential lifetime fitting there are two core algorithms within FLIMLib: The first is a triple integral method that does a very fast estimate of a single exponential lifetime component. The second is a Levenberg-Marquardt algorithm or LMA that uses an iterative, least-squares-minimization approach to generate a fit. This works with single, double and triple exponential models, as well as stretched exponential. There is also code to perform 'global' analysis over a number of signals symultaneously (e.g. over an image), where the lifetimes can be considered constant across the data set, but the amplitudes are allowed to vary for each signal. There is also a completely generic global analysis function. A third algorithm is available to perform phasor analysis.

In addition there is a non-negative linear least squares algorithm that is useful for spectral unmixing in combined spectral lifetime imaging (SLIM).

The code is written in C89 compatible C and is thread safe for fitting multiple pixels concurrently.  Several files are provided as wrappers to call this library from Java code: EcfWrapper.c and .h provide a subset of function calls used by FLIMJ.

For further details, see:
    https://flimlib.github.io/
    
If you are familiar with the program TRI2, that uses FLIMLib, this screenshot may help you to understand the meaning of the parameters.

@image html params_in_tri2.png
@image latex params_in_tri2.jpg "How some FLIMLib paramters are used in TRI2." width=15cm
    
    
# Library Contents:

Directory            | Contents
---------            | --------
src	                 | source files
src/main/c           | The source files for the  Curve library
src/flimlib-cmd/c | The source files for the stand alone executable wrapper for the library
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

    cp flimlib-cmd ../test_files
    
Run the program with the test files

    cd ../test_files
    ./flimlib-cmd test.ini transient.dat
    
    
\author Paul Barber
\copyright Creative Commons BY-SA 2014
    
*/
