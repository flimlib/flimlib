#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 2.6)

project(slim-curve-test)

include_directories("../slim-curve/src/main/c", ".")

# The current library
set(SLIM_CURVE_COMMON   ../slim-curve/src/main/c/EcfGlobal.c
                        ../slim-curve/src/main/c/EcfSingle.c
                        ../slim-curve/src/main/c/EcfSPA.c
                        ../slim-curve/src/main/c/EcfUtil.c
                        ../slim-curve/src/main/c/EcfWrapper.c
                        ../slim-curve/src/main/c/GCI_Lsqnonneg.c
                        ../slim-curve/src/main/c/GCI_Phasor.c
)

# Old library for comparison
set(SLIM_CURVE_OLD  	OldEcfGlobal.c
                        OldEcfSingle.c
)


# Additional files for the stand alone program
set(SLIM_CURVE_MAIN     main.c
                        json.c
                        json_helper.c
                        output.c
                        parser.c
                        tests.c
)

# Build a test program
add_executable(slim-curve-test ${SLIM_CURVE_COMMON} ${SLIM_CURVE_MAIN} ${SLIM_CURVE_OLD})

if(UNIX)
target_link_libraries(slim-curve-test m) #link the math library
endif(UNIX)
