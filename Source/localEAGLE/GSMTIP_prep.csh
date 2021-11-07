#!/bin/bash
rm *.o
#ifort -c -check all -g -debug full -traceback -132 -assume byterec mo_bas_gsm.f90
#ifort -c -check all -g -debug full -traceback -132 -assume byterecl *.for

#ifort -c -fpe0 -assume byterecl -132 *.for
ifort -c -O1 -assume byterecl -132 mo_gsm_const.f90
ifort -c -O1 -assume byterecl -132 mo_ham_gsm.f90
ifort -c -O1 -assume byterecl -132 *.for
ifort *.o
echo 'ifort end'
