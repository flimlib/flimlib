  $ cp "$TESTDIR/../target/build/bin/flimlib-cmd" "$TESTDIR"/../test_files/* .

  $ ./flimlib-cmd settings1.ini data.dat
  Settings: settings1.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.731979 X2 8.997931
  LMA return value 19
  LMA fitted A 11807.464844 T 2.268578 Z 252.651108 X2 7.323193

  $ ./flimlib-cmd settings2.ini data.dat
  Settings: settings2.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.731979 X2 8.997931
  LMA return value 21
  LMA fitted A1 8786.699219 T1 2.841778 A2 5456.799316 T2 0.779751 Z 100.473969 X2 1.210896

  $ ./flimlib-cmd settings3.ini data.dat
  Settings: settings3.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.731979 X2 8.997931
  LMA return value 34
  LMA fitted A1 7850.844238 T1 3.024728 A2 5444.664551 T2 1.022018 A3 4237.820312 T3 0.162697 Z 71.399940 X2 1.133069

  $ ./flimlib-cmd settings4.ini data.dat
  Settings: settings4.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.731979 X2 8.997931
  LMA return value 10
  LMA fitted A 16494.* T 2.6487* H 1.3779* Z -2194.* X2 204.4* (glob)

  $ ./flimlib-cmd settings5.ini data.dat
  Settings: settings5.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.731979 X2 8.997931
  Phasor return value 0
  Phasor fitted T 2.330856 Z 208.731979 X2 7.946539
