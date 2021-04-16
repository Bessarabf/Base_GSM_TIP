! 26/02/2020 approx co2 dencity

! at 80 km anCO2(1)=360 ppm
! at higher altitude cO2 density are calculated in diffusion
! equilibrium approximation.
! 
! anCO2 - CO2 density 1/cm-3
! an1   - O2  density 1/cm-3
! an2   - N2  density 1/cm-3
! an3   - O   density 1/cm-3
! an6   - temperature, K
! eddyco - 2d eddy diffusion coefficient

      subroutine co2den_eddy(anco2,an1,an2,an3,an6,eddyco,
     *                   rads,rp,g,n,n1,n2)
      dimension an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n),eddyco(n,n1)
     *         ,an6(n1,n2,n),anco2(n1,n2,n),rads(n),g(n),rp(n)
      data amco2/73.04e-24/
c
      do i=1,n1
        do j=1,n2
           anCO2(i,j,1)=360.e-6*an2(i,j,1)
        end do
      end do
      call bardif_eddy(anCO2,an1,an2,an3,an6,eddyco,rp,g,amCO2,
     *                 n,n1,n2,2)
      return
      end

! density approcsimation in diffusion equilibrium 
!  
! eddyco - 2d eddy dif coefficient
!
! an1   - O2  density 1/cm-3
! an2   - N2  density 1/cm-3
! an3   - O   density 1/cm-3
! an6   - temperature, K

      subroutine bardif_eddy(an,an1,an2,an3,an6,eddyco,rp,g,am,
     *                       n,n1,n2,l)
      USE mo_gsm_const, ONLY:amO2,amN2,amO,bk
      dimension an(n1,n2,n),an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n)
     *         ,an6(n1,n2,n),eddyco(n,n1),rp(n),g(n)
      data ves/0.5/
      f2=bk/am
      do 1 i=1,n1
       do 2 j=1,n2
        do 3 k=l,n
          tn=an6(i,j,k-1)
          tv=an6(i,j,k)
          hn=f2*tn/g(k-1)
          hv=f2*tv/g(k)
          sumV=(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))
          sumN=(an1(i,j,k-1)+an2(i,j,k-1)+an3(i,j,k-1))
          amsV=(an1(i,j,k)*amO2+an2(i,j,k)*amN2+an3(i,j,k)*amO)/
     *         sumV
          amsN=(an1(i,j,k-1)*amO2+an2(i,j,k-1)*amN2+an3(i,j,k-1)*amO)/
     *         sumN
          hsrV=tv*bk/(amsV*g(k))
          hsrN=tn*bk/(amsN*g(k-1))
!
          cMolv=3.e17/sumV*sqrt(tv)
          cMoln=3.e17/sumN*sqrt(tn)

          Hn=(cMoln*hn+eddyco(k-1,i)*hsrn)/(cMoln+eddyco(k-1,i))
          Hv=(cMolv*hv+eddyco(k,i)*hsrv)/(cMolv+eddyco(k,i))
          an(i,j,k)=an(i,j,k-1)*tn/tv*exp(-rp(k-1)*0.5*(1./Hn+1./Hv))

    3   continue
    2  continue
    1 continue
      return
      end

