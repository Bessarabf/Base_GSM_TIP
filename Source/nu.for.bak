c   nu:  lambet-tube-trubka-cyclt2
c   nu:  rik-lambet-tube-trubka-cyclt2
c   nu:  pitn-alga-tube-trubka-cyclt2
c   nu:  pqji-piqj- alga-tube-trubka-cyclt2
c   nu:  piqj- alga-tube-trubka-cyclt2
c
!      real function nu(ne,ti,tn)
!      data s1/.37/,s2/.38/,s3/-1.5/
!      goto(1,1,4),ne
!    1 continue
!        t=ti+tn
!        goto(2,3),ne
!    2   continue
!          a=t**s1
!          goto5
!    3   continue
!          a=t**s2
!          goto5
!    4   continue
!          a=ti**s3
!    5   continue
!        nu=a
!        return
!      end
! Updated Mihail Melnik from Polar Geophysical Institute 06.05.2021
      REAL FUNCTION NU(Ne,Ti,Tn)
      IMPLICIT NONE
      REAL  T
      REAL Ti , Tn  ! , intent(in)
      INTEGER Ne    ! , intent(in)
!      DATA s1/.37/ , s2/.38/ , s3/ - 1.5/
      T = Ti + Tn
      IF ( Ne==1 ) NU = T **(0.37)
      IF ( Ne==2 ) NU = T **(0.38)
!      IF ( Ne==3 ) NU = Ti**(-1.5)
      IF ( Ne==3 ) NU = 1.0 / sqrt ( Ti * Ti * Ti)
      END