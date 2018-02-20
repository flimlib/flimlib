  $ cp "$TESTDIR/../build/slim-curve-test" "$TESTDIR"/../test_files/* .

  $ (./slim-curve-test ./original.json) | grep '^[\s\S]*{' > ref.json
  
  $ ./slim-curve-test ./ref.json
  *
  
