c     . . . energy and flux of electron precipitations
      subroutine iacflo_AE(utn,sole,solu,bmpz,bmpy,mas,
     *                  csol,vsol,fa0,pkp0,ap0,ae0,dst0,al0,au0,
     *                  fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps,
     *                  AEj2)
      dimension sole(nse),solu(nsu),mas(10),solen(nse),ps(10)
! . . . soft electron (cusp) 70 deg
!     ps(1 )=1.e8
      ps(1 )=5.e8
!     ps(1 )=5.e9
!     ps(1 )=2.5e9
!     ps(1 )=1.5e9
      ps(2 )=1.
      ps(3 )=0.20e3
! . . . hard electron (auroral) 70 deg
!     ps(4 )=5.e8
      ps(4 )=1.e7  ! flux
      ps(5 )=1.0   ! param gamma
      ps(6 )=3.e3  ! characteristic energy
!     ps(6 )=5.e3
! . . .soft electrons 80 deg
      ps(7 )=0.e8  ! first zone flux
!     ps(7 )=5.e8
!     ps(7 )=4.e9
      ps(8 )=0.1e3 ! characteristic energy
!     ps(8 )=0.05e3
      ps(9 )=0.e8  ! second zone flux
!     ps(9 )=3.e8
      ps(10)=0.05e3 ! characteristic energy
! . . .geomagnetic indexes ...
      if(mas(2).eq.0) then
!       . . .from input data
        fa=fa0
        pkp=pkp0
        ap=ap0
        ae=ae0
        dst=dst0
        al=al0
        au=au0
!       . . .calculating
      else
        fa=1
        pkp=1
        ap=1
        ae=1
        dst=1
        al=1
        au=1
      end if
!!!!!!!!!!!     storm ..............
!    . . . «¿¬»—»ÃŒ—“‹ ¬€—€œ¿Õ»… Œ“ Kp
!1
!         ps(4)=ps(4)*(-3.5+6.5*pkpj2)
!1
!2
!         ps(4)=ps(4)*1.43*pkpj2
!2
!3
         ps(4)=ps(4)*(1.0+0.004*AE0)
!3
!         ps(6)=5.e3
         ps(6)=2.58e3+0.003e3*AE0
!!!!!!!!!!!     storm ..............
      return
      end
