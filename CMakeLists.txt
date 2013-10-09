#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 2.6)

project(slim-curve)

include_directories("src/main/c")

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
set(SLIM_CURVE_MAIN     src/main/c/main.c
)

# Aditional files for a library
set(SLIM_CURVE_LIB      src/main/c/loci_curvefitter_SLIMCurveFitter.c
                        src/main/c/loci_slim_SLIMCurve.c
)

# Build a dll
#add_library(slim-curve SHARED ${SLIM_CURVE_COMMON} ${SLIM_CURVE_LIB})

# Build a stand alone program
add_executable(slim-curve-cmd ${SLIM_CURVE_COMMON} ${SLIM_CURVE_MAIN})
target_link_libraries(slim-curve-cmd m) #link the math library
