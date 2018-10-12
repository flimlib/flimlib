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
sh travis-build.sh $encrypted_58cee4862e74_key $encrypted_58cee4862e74_iv

exit_code=$?

# NB: We skip tests when doing a release, because:
# A) The target folder structure differs, and cram hardcodes it; and
# B) The release already happened and was deployed by now. ;-)
if [ ! -f ./target/checkout/release.properties ]
then
  # Run the unit tests
  "$CRAM" ./tests

  exit_code=$((exit_code | $?))
fi

ls -lR ./target/

exit $exit_code
# TODO: Maven deploy artifacts
