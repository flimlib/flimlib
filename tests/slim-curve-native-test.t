  $ cp "$TESTDIR/../build/slim-curve-test" "$TESTDIR"/../test_files/* .

# greps changed json string from output as new reference
  $ (./slim-curve-test ./original.json) | grep '^[\s\S]*{' > ref.json

  # finds FAILed tests from output
  $ (./slim-curve-test ./ref.json) | perl -0777 -ne "while(m/\nTEST((.*\S.*)\n)+FAIL((.*\S.*)\n)+/g){print $&}"

