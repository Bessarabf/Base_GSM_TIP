:: ifort /c /assume:byterecl /debug:full mo_bas_gsm.f90
:: ifort /c /assume:byterecl /Od /check:all /debug:full *.for 

ifort /c Source\mo_bas_gsm.f90
::ifort /c /assume:byterecl /Od /extend_source:80 /fpp /Qopenmp  /traceback /check:bounds  *.for 
ifort /c /assume:byterecl /extend_source:80 /fpp /traceback /check:bounds  Source\*.for 
:: ifort *.obj /link /OUT:"GSM_TIP.exe" 
ifort *.obj
del *.obj
del mo_bas_gsm.mod
aepp.exe