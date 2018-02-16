#!/bin/sh

# Install needed tools
if which cram; then
  CRAM=cram
else
  pip install --user cram
  CRAM=/cram
  CRAM=$HOME/.local/bin/cram
fi

# Prepare the build environment
rm -rf build
mkdir -p build
cd build
cmake ..

# Build the code
make

# Run the unit tests
"$CRAM" ../tests

# TODO: Maven deploy artifacts
