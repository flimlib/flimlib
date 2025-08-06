#!/bin/sh

# Deploy script with aggregated multi-platform native artifacts.
# Runs on the Linux deploy node after downloading artifacts from all platforms.

echo '--> Moving downloaded platform artifacts into place:'
mkdir -p target

# Drop in the non-matching-platform native JARs.
for jar in \
  */flimlib-*-natives-osx_*.jar \
  */flimlib-*-natives-windows_*.jar
do
  if [ -f "$jar" ]; then mv -fv "$jar" target/; fi
done

echo
echo 'Available artifacts:'
ls -la target/
echo

# Run a full build to generate SWIG classes needed for javadoc,
# then deploy all JARs (including the downloaded platform artifacts).
curl -fsLO https://raw.githubusercontent.com/scijava/scijava-scripts/main/ci-build.sh
BUILD_ARGS='-Pattach-prebuilt-artifacts' sh ci-build.sh
