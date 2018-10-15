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

# Deploy artifacts
# Lifted from https://github.com/imagej/imagej-launcher/blob/f14435e80acbe7c84d52695a4794afb570ab65c8/.travis/build.sh
# Get GAV
groupId="$(sed -n 's/^  <groupId>\(.*\)<\/groupId>$/\1/p' pom.xml)"
groupIdForURL="$(echo $groupId | sed -e 's/\./\//g')"
artifactId="$(sed -n 's/^ <artifactId>\(.*\)<\/artifactId>$/\1/p' pom.xml)"
version="$(sed -n 's/^  <version>\(.*\)<\/version>$/\1/p' pom.xml)"

artifactPaths="target/checkout/target/*.jar"

if [ "$TRAVIS_OS_NAME" = "linux" ]
then
  classifier="linux_64"
else
  classifier="osx_64"
fi

if [ "$TRAVIS_SECURE_ENV_VARS" = true \
  -a "$TRAVIS_PULL_REQUEST" = false \
  -a -f "target/checkout/release.properties" ]
then
  echo "== Deploying binaries =="

  # Check if a release has been deployed for that version
  folderStatus=$(curl -s -o /dev/null -I -w '%{http_code}' http://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/)
  if [ "$folderStatus" != "200" ]
  then
    exit $exit_code
  fi

  for artifactPath in artifactPaths; do
    artifactFileName="${artifactPath##*/}"
    extension="${artifactPath##*.}"
    # Check if the launcher for that version has already been deployed
    fileStatus=$(curl -s -o /dev/null -I -w '%{http_code}' http://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/$artifactFileName)
    if [ "$fileStatus" != "200" ]
    then
      mvn deploy:deploy-file\
        -Dfile="$artifactPath"\
        -DrepositoryId="imagej.releases"\
        -Durl="dav:https://maven.imagej.net/content/repositories/releases"\
        -DgeneratePom="false"\
        -DgroupId="$groupId"\
        -DartifactId="$artifactId"\
        -Dversion="$version"\
        -Dclassifier="$classifier"\
        -Dpackaging="$extension"
    fi
  done
fi

exit $exit_code
