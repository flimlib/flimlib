#
# CMakeLists.txt
#

cmake_minimum_required(VERSION 2.6)

project(flimlib)

# CXX Build options
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MACOSX_RPATH 1)

# Common output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories("src/main/c" "src/main/cpp" "src/flimlib-cmd/c/iniparser/src" "src/test/c")

# Common files we always need
set(FLIMLIB_COMMON   src/main/c/EcfGlobal.c
                     src/main/c/EcfSingle.c
                     src/main/c/EcfSPA.c
                     src/main/c/EcfUtil.c
                     src/main/c/EcfWrapper.c
                     src/main/c/GCI_Lsqnonneg.c
                     src/main/c/GCI_Phasor.c
)

# Additional files for the stand alone program
set(FLIMLIB_MAIN     src/flimlib-cmd/c/main.c
                     src/flimlib-cmd/c/iniparser/src/iniparser.c
                     src/flimlib-cmd/c/iniparser/src/dictionary.c
)

# Additional files for the stand alone testing
set(FLIMLIB_TEST     src/test/c/main.c
                     src/test/c/json.c
                     src/test/c/json_helper.c
                     src/test/c/output.c
                     src/test/c/parser.c
                     src/test/c/tests.c
)

# Aditional files for a library
set(FLIMLIB_LIB      src/main/c/flimlib.def
)

# Aditional files for a stand alone program based on the C++ interface
set(FLIMLIB_CPP      src/flimlib-cmd/cpp/main.cpp
                     src/main/cpp/flimlib.hpp
                     src/flimlib-cmd/c/iniparser/src/iniparser.c
                     src/flimlib-cmd/c/iniparser/src/dictionary.c
)

# Build a dll
add_library(flimlib SHARED ${FLIMLIB_COMMON} ${FLIMLIB_LIB})

# Build a static lib
add_library(flimlib-static STATIC ${FLIMLIB_COMMON})

# Build a dll test program
add_executable(flimlib-dll-test  ${FLIMLIB_MAIN})
target_link_libraries(flimlib-dll-test flimlib)

# Build a stand alone program
add_executable(flimlib-cmd ${FLIMLIB_COMMON} ${FLIMLIB_MAIN})

# Build a stand alone c++ program
add_executable(flimlib-cmd-cpp ${FLIMLIB_COMMON} ${FLIMLIB_CPP})

# Build a test program
add_executable(flimlib-test ${FLIMLIB_COMMON} ${FLIMLIB_TEST})

if(UNIX)
  # Link the math library
  target_link_libraries(flimlib m)
  target_link_libraries(flimlib-cmd m)
  target_link_libraries(flimlib-test m)
endif(UNIX)

# Setup swig to build a java lib
set(CMAKE_SWIG_FLAGS -package flimlib)
# simply to avoid warnings, see https://github.com/CCPPETMR/SIRF/issues/242
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13.5") 
    # policy introduced in CMake 3.13
    cmake_policy(SET CMP0086 NEW)
    cmake_policy(SET CMP0078 NEW)
endif()

find_package(SWIG REQUIRED)
find_package(Java REQUIRED)
find_package(JNI REQUIRED)
include(UseJava)
include(UseSWIG)
include_directories(${JNI_INCLUDE_DIRS})

# direcotries where the generated java classes and c wrapper are kept
set( CMAKE_SWIG_OUTDIR ${PROJECT_SOURCE_DIR}/target/generated-sources/main/java/flimlib )
file(MAKE_DIRECTORY ${CMAKE_SWIG_OUTDIR})
set( SWIG_OUTFILE_DIR ${PROJECT_SOURCE_DIR}/target/generated-sources/main/cpp )
file(MAKE_DIRECTORY ${SWIG_OUTFILE_DIR})

# Build the C code into a dynamic library
set_property(SOURCE src/main/swig/FLIMLib.i PROPERTY CPLUSPLUS ON)
swig_add_library(
    flimlib-jni 
    TYPE SHARED
    LANGUAGE java
    SOURCES 
    	src/main/swig/FLIMLib.i
)
swig_link_libraries(
	flimlib-jni
	flimlib
)

# Detect target bitness and architecture
set( TARGET_NATIVE ${PROJECT_SOURCE_DIR}/target/natives )
math(EXPR BITS "8*${CMAKE_SIZEOF_VOID_P}")
set( OUT_DIR_W_BITS ${TARGET_NATIVE} ) #/x${BITS} )
file(MAKE_DIRECTORY ${OUT_DIR_W_BITS})

if(WIN32)
    set( TARGET_LIB_PREFIX "" )
else()
    set( TARGET_LIB_PREFIX "lib" )
endif()
# RUNTIME_ for dll output and LIBRARY_ for .so/.dylib/.jnilib output
# _DEBUG _RELEASE are for MSVC in particular
set_target_properties(
    flimlib
    PROPERTIES
        PREFIX                              "${TARGET_LIB_PREFIX}"
        RUNTIME_OUTPUT_DIRECTORY            "${OUT_DIR_W_BITS}"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG      "${OUT_DIR_W_BITS}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE    "${OUT_DIR_W_BITS}"
        LIBRARY_OUTPUT_DIRECTORY            "${OUT_DIR_W_BITS}"
)
set_target_properties(
    flimlib-jni
    PROPERTIES
        PREFIX                              "${TARGET_LIB_PREFIX}"
        RUNTIME_OUTPUT_DIRECTORY            "${OUT_DIR_W_BITS}"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG      "${OUT_DIR_W_BITS}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE    "${OUT_DIR_W_BITS}"
        LIBRARY_OUTPUT_DIRECTORY            "${OUT_DIR_W_BITS}"
)

######################debug utility: print all properties######################
# # Get all propreties that cmake supports (for debug purposes)
# execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)

# # Convert command output into a CMake list
# STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
# STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

# function(print_properties)
#     message ("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
# endfunction(print_properties)

# function(print_target_properties tgt)
#     if(NOT TARGET ${tgt})
#       message("There is no target named '${tgt}'")
#       return()
#     endif()

#     foreach (prop ${CMAKE_PROPERTY_LIST})
#         string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
#     # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
#     if(prop STREQUAL "LOCATION" OR prop MATCHES "^LOCATION_" OR prop MATCHES "_LOCATION$")
#         continue()
#     endif()
#         # message ("Checking ${prop}")
#         get_property(propval TARGET ${tgt} PROPERTY ${prop} SET)
#         if (propval)
#             get_target_property(propval ${tgt} ${prop})
#             message ("${tgt} ${prop} = ${propval}")
#         endif()
#     endforeach(prop)
# endfunction(print_target_properties)

# print_target_properties(flimlib-jni)
