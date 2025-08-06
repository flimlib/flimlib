#!/bin/sh

# Deploy script with aggregated multi-platform native artifacts.
# Runs on the deploy node after downloading artifacts from all platforms.

echo "--> Moving downloaded artifacts into place:"
mkdir -p target
for jar in */*.jar; do
  mv -fv "$jar" target/
done

echo
echo "Available artifacts:"
ls -la target/
echo

# Run the build in multiplatform deploy mode.
# This skips compilation etc. and deploys all JARs.
curl -fsLO https://raw.githubusercontent.com/scijava/scijava-scripts/main/ci-build.sh
sh ci-build.sh -Pdeploy-multiplatform
