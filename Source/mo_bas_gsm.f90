MODULE mo_bas_gsm
  
  INTEGER, PARAMETER :: ids0 = 72    ! number of longitudes
  INTEGER, PARAMETER :: its0 = 37    ! number of latitudes
  INTEGER, PARAMETER :: nh0 = 30     ! number of GSMTIP levels
  INTEGER, PARAMETER :: mlev =10     ! number of hamtoGSMTIP levels
! massives
  DIMENSION qdis(2,its0,ids0,nh0)
! constant section !!!!!!!!!!!!!!!!!!!
  real,parameter :: pi=3.1415926 
  real,parameter :: om=7.2722e-5
  real,parameter :: g0=980.665
  real,parameter :: re=6371.02e5 
  real,parameter :: bk=1.38064e-16 
  real,parameter :: amo2=53.12e-24,amn2=46.51e-24,amo=26.56e-24
  real,parameter :: amno=49.82e-24,amn=23.26e-24
END MODULE mo_bas_gsm