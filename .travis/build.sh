#!/bin/sh

# Install needed tools
if which cram; then
  CRAM=cram
else
  # Need to install cram
  if ! which pip; then
    # Need to install pip
    if [ $(uname) = "Darwin" ]; then
      brew install python
    else
      echo "No cram, no pip. Cannot continue."
      exit 1
    fi
  fi
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
