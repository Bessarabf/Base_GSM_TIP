!  Nusinov model UV
      subroutine nus_uv(solu,nsu)
      PARAMETER(nsu05=12)	! half of nsu
      dimension B0(nsu05-1),B1(nsu05-1),alyam(nsu05-1)
      dimension solu(2*nsu)

!  Nusinov et al., A Model of Fluxes of Solar Ultraviolet Irradiance
!  Geom. and aeronom,2019, 59, N 3, pp 284-290.
!
!  Lam1     Lam2          B0          B1   
! 1220.0   1250.0      0.027404       0.015021
! 1250.0   1270.0      0.009701       0.008101
! 1270.0   1310.0      0.061811       0.019138
! 1310.0   1350.0      0.031552       0.034689
! 1350.0   1380.0      0.040150       0.006783
! 1380.0   1500.0      0.282535       0.055582
! 1500.0   1550.0      0.259750       0.051551
! 1550.0   1630.0      0.840660       0.087146
! 1630.0   1670.0      0.806700       0.083010
! 1670.0   1720.0      1.752800       0.131150
! 1720.0   1760.0      2.392400       0.141150

      data B0 / 0.027404,0.009701,0.061811,0.031552,0.040150,
     *          0.282535,0.259750,0.840660,0.806700,1.752800,2.392400/
      data B1 / 0.015021,0.008101,0.019138,0.034689,0.006783,
     *          0.055582,0.051551,0.087146,0.083010,0.131150,0.141150/
!  wave length in nm (center of interval)
      data alyam /123.5,126,129,133,136.5,144,152.5,159,165,169.5,174/

      p_La=solu(1)*1.e-2 ! normalisation La for Nusinov model
      
      do k=2,nsu05
         solu(k) = (B0(k-1) + p_la*B1(k-1))*1.0e2 ! *10**9 phot/cm2/s
         solu(k+nsu05) = solu(k)*1.98648/alyam(k-1) ! erg/cm2/s
      end do
      return
      end
         
