#!/bin/sh

# Prepare the build environment
rm -rf build
mkdir -p build
cd build
cmake ..

# Build the code
make

# TODO: Maven deploy artifacts
