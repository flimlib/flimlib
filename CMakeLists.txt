#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 2.6)

project(slim-curve)

include_directories("src/main/c" "src/main/cpp" "src/slim-curve-cmd/c/iniparser/src")

# Common files we always need
set(SLIM_CURVE_COMMON   src/main/c/EcfGlobal.c
                        src/main/c/EcfSingle.c
                        src/main/c/EcfSPA.c
                        src/main/c/EcfUtil.c
                        src/main/c/EcfWrapper.c
                        src/main/c/GCI_Lsqnonneg.c
                        src/main/c/GCI_Phasor.c
)

# Additional files for the stand alone program
set(SLIM_CURVE_MAIN     src/slim-curve-cmd/c/main.c
                        src/slim-curve-cmd/c/iniparser/src/iniparser.c
                        src/slim-curve-cmd/c/iniparser/src/dictionary.c
)

# Aditional files for a library
set(SLIM_CURVE_LIB      #src/main/c/loci_curvefitter_SLIMCurveFitter.c
                        #src/main/c/loci_slim_SLIMCurve.c
                        src/main/c/slim-curve.def
)

# Aditional files for a stand alone program based on the C++ interface
set(SLIM_CURVE_CPP      src/slim-curve-cmd/cpp/main.cpp
                        src/main/cpp/slim-curve.hpp
                        src/slim-curve-cmd/c/iniparser/src/iniparser.c
                        src/slim-curve-cmd/c/iniparser/src/dictionary.c
)

# Build a dll
add_library(slim-curve SHARED ${SLIM_CURVE_COMMON} ${SLIM_CURVE_LIB})
target_link_libraries(slim-curve m)

# Build a dll test program
add_executable(slim-curve-dll-test  ${SLIM_CURVE_MAIN})
target_link_libraries(slim-curve-dll-test slim-curve)

# Build a stand alone program
add_executable(slim-curve-cmd ${SLIM_CURVE_COMMON} ${SLIM_CURVE_MAIN})

# Build a stand alone c++ program
add_executable(slim-curve-cmd-cpp ${SLIM_CURVE_COMMON} ${SLIM_CURVE_CPP})

if(UNIX)
target_link_libraries(slim-curve-cmd m) #link the math library
endif(UNIX)
