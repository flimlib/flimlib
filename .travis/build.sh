#!/bin/sh

# Install needed tools
if which brew; then
  brew install swig
fi
if which cram; then
  CRAM=cram
else
  # Need to install cram
  if which pip; then
    PIP=pip
  elif which pip2; then
    PIP=pip2
  else
    echo "No cram, no pip. Cannot continue."
    exit 1
  fi
  "$PIP" install --user cram
  if which cram; then
    CRAM=cram
  else
    for dir in \
      /usr/local/bin \
      "$HOME/Library/Python/2.7/bin" \
      "$HOME/.local/bin"
    do
      test -f "$dir/cram" && CRAM="$dir/cram"
    done
  fi
  if [ -z "$CRAM" ]; then
    echo "Cram purportedly installed, but cannot find it."
    "$PIP" show -f cram
    exit 2
  fi
fi

curl -fsLO https://raw.githubusercontent.com/scijava/scijava-scripts/master/travis-build.sh
sh travis-build.sh

exit_code=$?

# Run the unit tests
"$CRAM" ./tests

exit_code=$((exit_code | $?))

ls -l ./target/
ls -l ./target/natives/

exit $exit_code
# TODO: Maven deploy artifacts
