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
