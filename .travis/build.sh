#!/bin/bash

# Install needed tools
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

# ls -lR ./target/

# Deploy artifacts
# Lifted from https://github.com/imagej/imagej-launcher/blob/f14435e80acbe7c84d52695a4794afb570ab65c8/.travis/build.sh
extraArtifactPaths="./target/checkout/target/*.jar"

if [ "$TRAVIS_OS_NAME" = "linux" ]
then
  classifier="natives-linux_64"
else
  classifier="natives-osx_64"
fi

if [ "$TRAVIS_SECURE_ENV_VARS" = true \
  -a "$TRAVIS_PULL_REQUEST" = false \
  -a -f "target/checkout/release.properties" ]
then
  echo "== Deploying binaries =="

  # Get GAV
  gav=$(mvn -q \
    -Dexec.executable=echo \
    -Dexec.args='${project.groupId}:${project.artifactId}:${project.version}' \
    --non-recursive \
    exec:exec)
  groupId="$(echo $gav | cut -d ":" -f 1)"
  groupIdForURL="$(echo $groupId | sed -e 's/\./\//g')"
  artifactId="$(echo $gav | cut -d ":" -f 2)"
  version="$(echo $gav | cut -d ":" -f 3)"

  # Check if a release has been deployed for that version
  folderStatus=$(curl -s -o /dev/null -I -w '%{http_code}' https://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/)
  if [ "$folderStatus" != "200" ]
  then
    exit $exit_code
  fi

  for artifactPath in $extraArtifactPaths; do
    fileName="${artifactPath##*/}"
    # Skip the non-classified artifacts
    if [[ ! "$fileName" =~ "$classifier" ]]
    then
      continue
    fi
    extension="${fileName##*.}"
    # Check if the launcher for that version has already been deployed
    fileStatus=$(curl -s -o /dev/null -I -w '%{http_code}' https://maven.imagej.net/content/repositories/releases/$groupIdForURL/$artifactId/$version/$fileName)
    if [ "$fileStatus" != "200" ]
    then
      if [ ! -z "$mainFile" ]
      then
        files="$files,$mainFile"
        types="$types,$mainType"
        classifiers="$classifiers,$mainClassifier"
      fi
      mainFile="$artifactPath"
      mainType="$extension"
      mainClassifier="$classifier"
    fi
  done
  if [ ! -z "$mainFile" ]
  then
    deployCmd="mvn deploy:deploy-file\
      -Dfile=\"$mainFile\"\
      -DrepositoryId=\"imagej.releases\"\
      -Durl=\"dav:https://maven.imagej.net/content/repositories/releases\"\
      -DgeneratePom=\"false\"\
      -DgroupId=\"$groupId\"\
      -DartifactId=\"$artifactId\"\
      -Dversion=\"$version\"\
      -Dclassifier=\"$mainClassifier\"\
      -Dpackaging=\"$mainType\""
    if [ ! -z "$files" ]
    then
      deployCmd="$deployCmd -Dfiles=\"$files\" -Dclassifiers=\"$classifiers\" -Dtypes=\"$types\""
    fi
    eval "$deployCmd"
  fi
  exit_code=$((exit_code | $?))
fi

exit $exit_code
