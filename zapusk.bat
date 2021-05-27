:: ifort /c /assume:byterecl /debug:full mo_bas_gsm.f90
:: ifort /c /assume:byterecl /Od /check:all /debug:full *.for 

ifort /c Source\mo_bas_gsm.f90
:: ifort /c /assume:byterecl /Od /extend_source:80 /fpp /traceback /check:bounds  /check:all Source\*.for 
:: ifort /c /assume:byterecl /extend_source:80 /fpp /traceback  /Qfp-speculation=off /fpe=0      Source\*.for
   ifort /c /assume:byterecl /extend_source:80 /fpp /traceback  /Qfp-speculation=off /fpe=0 /Od  Source\*.for 
ifort *.obj /link /OUT:"GSM_TIP.exe" 
del *.obj
del mo_bas_gsm.mod
GSM_TIP.exe >output_d.txt
::GSM_TIP.exe