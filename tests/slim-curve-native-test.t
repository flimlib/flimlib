  $ cp "$TESTDIR/../build/slim-curve-test" "$TESTDIR"/../test_files/* .

# greps changed json string from output as new reference
  $ (./slim-curve-test ./original.json) | grep '^[\s\S]*{' > ref.json

# This test should pass
  $ ./slim-curve-test temp.json | grep -Pzo "\n*TEST(.+\n)*FAIL(.+\n)*"; test $? -eq 1

# greps FAILed tests from output and reverse return code of grep
  $ ./slim-curve-test ref.json | grep -Pzo "\n*TEST(.+\n)*FAIL(.+\n)*"; test $? -eq 1
