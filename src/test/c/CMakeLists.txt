#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 2.6)

project(flimlib-test)

include_directories("../flimlib/src/main/c", ".")

# The current library
set(FLIMLIB_COMMON   ../flimlib/src/main/c/EcfGlobal.c
                     ../flimlib/src/main/c/EcfMultiple.c
                     ../flimlib/src/main/c/EcfSPA.c
                     ../flimlib/src/main/c/EcfSingle.c
                     ../flimlib/src/main/c/EcfUtil.c
                     ../flimlib/src/main/c/EcfWrapper.c
                     ../flimlib/src/main/c/GCI_Lsqnonneg.c
                     ../flimlib/src/main/c/GCI_Phasor.c
)


# Additional files for the stand alone program
set(FLIMLIB_MAIN     main.c
                     json.c
                     json_helper.c
                     output.c
                     parser.c
                     tests.c
)

# Build a test program
add_executable(flimlib-test ${FLIMLIB_COMMON} ${FLIMLIB_MAIN})

if(UNIX)
target_link_libraries(flimlib-test m) #link the math library
endif(UNIX)
