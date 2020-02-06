!  diffusion equilibrium 
!  
!
      subroutine bardif(an,an1,an2,an3,an6,ctd,rp,g,am,n,n1,n2,l)
      USE mo_bas_gsm, ONLY:amO2,amN2,amO,bk
      dimension an(n1,n2,n),an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n)
     *         ,an6(n1,n2,n),ctd(n),rp(n),g(n)
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

          Hn=(cMoln*hn+ctd(k-1)*hsrn)/(cMoln+ctd(k-1))
          Hv=(cMolv*hv+ctd(k)*hsrv)/(cMolv+ctd(k))
          an(i,j,k)=an(i,j,k-1)*tn/tv*exp(-rp(k-1)*0.5*(1./Hn+1./Hv))

    3   continue
    2  continue
    1 continue
      return
      end
