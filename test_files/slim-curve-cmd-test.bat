@echo off
echo Test 1 - Mono
call slim-curve-cmd settings1.ini data.dat > console_output.txt
findstr /C:"RLD estimate A 11515.900391 T 2.392674 Z 208.732056 X2 8.997932" console_output.txt || goto :fail
findstr /C:"LMA fitted A 11807.464844 T 2.268578 Z 252.651108 X2 7.323208" console_output.txt || goto :fail
echo Test 2 - Bi
call slim-curve-cmd settings2.ini data.dat > console_output.txt
findstr /C:"LMA fitted A1 8786.449219 T1 2.841830 A2 5456.925293 T2 0.779802 Z 100.463684 X2 1.210934" console_output.txt || goto :fail
echo Test 3 - Tri
call slim-curve-cmd settings3.ini data.dat > console_output.txt
findstr /C:"LMA fitted A1 7851.750488 T1 3.024554 A2 5444.777344 T2 1.021768 A3 4244.931641 T3 0.162466 Z 71.422325 X2 1.133109" console_output.txt || goto :fail
echo Test 4 - Stretched
call slim-curve-cmd settings4.ini data.dat > console_output.txt
findstr /C:"LMA fitted A 16494.882813 T 2.648694 H 1.378027 Z -2194.759277 X2 204.380937" console_output.txt || goto :fail
echo Test 5 - Phasor
call slim-curve-cmd settings5.ini data.dat > console_output.txt
findstr /C:"Phasor fitted T 2.330857 Z 208.732056 X2 7.946540" console_output.txt || goto :fail
echo ****** All tests passed ******
pause
goto :eof

:fail
echo Failed!
pause
goto :eof