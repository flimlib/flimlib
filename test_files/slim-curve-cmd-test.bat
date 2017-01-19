@echo off
if not exist slim-curve-cmd.exe goto :file_error
echo Test 1 - Mono
call slim-curve-cmd settings1.ini data.dat > console_output.txt
findstr /R /C:"RLD estimate A 11515.900[0-9]* T 2.392[0-9]* Z 208.731[0-9]* X2 8.997[0-9]*" console_output.txt || goto :fail
findstr /R /C:"LMA fitted A 11807.464[0-9]* T 2.268[0-9]* Z 252.651[0-9]* X2 7.323[0-9]*" console_output.txt || goto :fail
echo Test 2 - Bi
call slim-curve-cmd settings2.ini data.dat > console_output.txt
findstr /R /C:"LMA fitted A1 8786.699[0-9]* T1 2.841[0-9]* A2 5456.799[0-9]* T2 0.779[0-9]* Z 100.473[0-9]* X2 1.210[0-9]*" console_output.txt || goto :fail
echo Test 3 - Tri
call slim-curve-cmd settings3.ini data.dat > console_output.txt
findstr /R /C:"LMA fitted A1 7850.844[0-9]* T1 3.024[0-9]* A2 5444.664[0-9]* T2 1.022[0-9]* A3 4237.820[0-9]* T3 0.162[0-9]* Z 71.399[0-9]* X2 1.133[0-9]*" console_output.txt || goto :fail
echo Test 4 - Stretched
call slim-curve-cmd settings4.ini data.dat > console_output.txt
findstr /R /C:"LMA fitted A 16494.630[0-9]* T 2.648[0-9]* H 1.377[0-9]* Z -2194.706[0-9]* X2 204.412[0-9]*" console_output.txt || goto :fail
echo Test 5 - Phasor
call slim-curve-cmd settings5.ini data.dat > console_output.txt
findstr /R /C:"Phasor fitted T 2.330[0-9]* Z 208.731[0-9]* X2 7.946[0-9]*" console_output.txt || goto :fail
echo ****** All tests passed ******
pause
goto :eof

:file_error
echo Could not find slim-curve-cmd-exe, place it in this folder.
pause
goto :eof

:fail
echo Output was...
type console_output.txt
echo That FAILED!
pause
goto :eof