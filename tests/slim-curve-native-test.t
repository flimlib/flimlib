  $ cp "$TESTDIR/../build/slim-curve-test" "$TESTDIR"/../test_files/* .

  $ (./slim-curve-test ./original.json) | grep '^[\s\S]*{' > ref.json

  $ cat ref.json
  * (glob)

  $ ./slim-curve-test ./ref.json
  *, 0 FAILURES (glob)
