  $ cp "$TESTDIR/../target/build/bin/flimlib-cmd" "$TESTDIR"/../test_files/* .

  $ ./flimlib-cmd settings1.ini data.dat
  Settings: settings1.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.731* X2 8.9979* (glob)
  LMA return value [12][0-9] (re)
  LMA fitted A 11807.4* T 2.2685* Z 252.6* X2 7.3231* (glob)

  $ ./flimlib-cmd settings2.ini data.dat
  Settings: settings2.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.7319* X2 8.99793* (glob)
  LMA return value [12][0-9] (re)
  LMA fitted A1 8786.69* T1 2.841778 A2 5456.799* T2 0.77975* Z 100.473* X2 1.210* (glob)

  $ ./flimlib-cmd settings3.ini data.dat
  Settings: settings3.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.7319* X2 8.99793* (glob)
  LMA return value 3? (glob)
  LMA fitted A1 785?.* T1 3.02* A2 544?.* T2 1.02* A3 42??.* T3 0.162* Z 71.* X2 1.133* (glob)

  $ ./flimlib-cmd settings4.ini data.dat
  Settings: settings4.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.7319* X2 8.99793* (glob)
  LMA return value [1-5][0-9] (re)
  LMA fitted A 15745.* T 1.6161* H 1.3704* Z 0.0* X2 1.145* (glob)

  $ ./flimlib-cmd settings5.ini data.dat
  Settings: settings5.ini
  Transient Data: data.dat
  TRI2 TR Settings File
  RLD return_value 2
  RLD estimate A 11515.900391 T 2.392674 Z 208.7319* X2 8.99793* (glob)
  Phasor return value 0
  Phasor fitted T 2.330856 Z 208.7319* X2 7.9465* (glob)
