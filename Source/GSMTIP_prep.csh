#!/bin/bash
rm *.o
#ifort -c -check all -g -debug full -traceback -132 -assume byterec mo_bas_gsm.f90
#ifort -c -check all -g -debug full -traceback -132 -assume byterecl *.for

#ifort -c -O0 -assume byterecl -132 mo_bas_gsm.f90
#ifort -c -O0 -assume byterecl -132 *.for

ifort -c -assume byterecl -extend-source 80 -O1 -fp-speculation=off -fpe=0 mo_bas_gsm.f90
ifort -c -assume byterecl -extend-source 80 -O1 -fp-speculation=off -fpe=0 *.for
ifort *.o

echo 'ifort end'
