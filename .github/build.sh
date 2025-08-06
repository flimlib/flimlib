#!/bin/sh
curl -fsLO https://raw.githubusercontent.com/scijava/scijava-scripts/main/ci-build.sh
NO_DEPLOY=1 sh ci-build.sh

exit_code=$?
if [ $exit_code -ne 0 ]; then exit $exit_code; fi

# Run the unit tests.
if [ "$BUILD_OS" = Windows ]; then
  echo "[WARNING] Skipping prysk tests on Windows platform."
else
  python -m pip install --upgrade pip
  pip install prysk==0.12.2
  python -m prysk tests
fi
