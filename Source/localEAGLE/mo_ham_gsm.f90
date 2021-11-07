MODULE mo_ham_gsm
  
  INTEGER, PARAMETER :: ids0 = 72    ! number of longitudes
  INTEGER, PARAMETER :: its0 = 37    ! number of latitudes
  INTEGER, PARAMETER :: nh0 = 30    ! number of GSMTIP levels
  INTEGER, PARAMETER :: mlev =1     ! 10 number of hamtoGSMTIP levels

! dissociation massive
  DIMENSION qdis(2,its0,ids0,nh0)
! revert mass for HAM data:
  ! tn and density  
  DIMENSION gsmHAM(nh0,its0,ids0),dgsmHAM(nh0,its0,ids0)
  ! velocity
  DIMENSION  UgsmHam(nh0,its0,ids0),VgsmHam(nh0,its0,ids0)
  ! n(o), n(no), n(n)
  DIMENSION  dOgsm(nh0,its0,ids0),dNOgsm(nh0,its0,ids0),dNgsm(nh0,its0,ids0)
! particle ionization 
  DIMENSION zpiongsm(nh0,its0,ids0)
! JoulHeating 
  DIMENSION qJGSM(its0,ids0,nh0)

! Ion drag
  DIMENSION dragViGSM(its0,ids0,nh0),dragVjGSM(its0,ids0,nh0)

END MODULE mo_ham_gsm

