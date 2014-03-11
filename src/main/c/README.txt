/**
@mainpage SLIM-curve package for exponential curve fitting of spectral lifetime data.

SLIM Curve is an exponential curve fitting library used for Fluorescent
Lifetime Imaging or FLIM and Spectral Lifetime Imaging or SLIM.  It is based on
code developed by Paul Barber and his group at the Gray Institute for Radiation
Oncology & Biology and used for FLIM functionality in his TRI2 (Time Resolved
Imaging) software.  It is also used in the LOCI SLIM Plugin project.

There are two algorithms used for curve fitting within SLIM Curve: The first is
a triple integral method that does a very fast estimate of a single exponential
lifetime component. The second is a Levenberg-Marquardt algorithm or LMA that
uses an iterative, least-squares-minimization approach to generate a fit. This
works with single, double and triple exponential models, as well as stretched
exponential.

The code is written in C89 compatible C and is thread safe for fitting multiple
pixels concurrently.  Several files are provided as wrappers to call this
library from Java code:  EcfWrapper.c and .h provide a subset of function calls
used by SLIM Plugin, these may be invoked directly from Java using JNA.  In
addition there is a Java CurveFitter project that provides a wrapper to the
SLIM Curve code.  This invokes the C code using JNI, with
loci_curvefitter_SLIMCurveFitter.c and .h.

For further details, see:
    http://loci.wisc.edu/software/slim-curve

Directory Contents:<br>
src						source files<br>
src/main/c				The source files for the SLIM Curve library<br>
src/slim-curve-cmd/c	The source files for the stand alone executable wrapper for the library<br>
test_files				dat and ini settings file for testing<br>
src/main/c/doc			API documentation (Doxygen output)<br>

To build the stand alone program using CMake and gcc under Linux:<br>
Create a build folder, and cd to it<br>
&nbsp;&nbsp;&nbsp;&nbsp;mkdir build<br>
&nbsp;&nbsp;&nbsp;&nbsp;cd build<br>
Run CMake<br>
&nbsp;&nbsp;&nbsp;&nbsp;cmake ../CMakeLists.txt<br>
Run make<br>
&nbsp;&nbsp;&nbsp;&nbsp;make<br>
	
To Run the stand alone executable<br>
Copy the executable to the test_files folder for convenience<br>
&nbsp;&nbsp;&nbsp;&nbsp;cp slim-curve-cmd ../test_files<br>
Run the program with the test files<br>
&nbsp;&nbsp;&nbsp;&nbsp;cd ../test_files<br>
&nbsp;&nbsp;&nbsp;&nbsp;./slim-curve-cmd test.ini transient.dat<br>
*/
