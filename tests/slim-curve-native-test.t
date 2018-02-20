  $ cp "$TESTDIR/../build/slim-curve-test" "$TESTDIR"/../test_files/* .

  $ sed -n 's/^[\s\S]*{/{/w ref.json' <<< $(./slim-curve-test ./original.json)

  $ ./slim-curve-test ./ref.json
  