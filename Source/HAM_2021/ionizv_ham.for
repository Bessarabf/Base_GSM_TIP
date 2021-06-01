      subroutine ionizv_ham(i0,gam,e0,emin,emax,h,par,parj,peion,
     *                      kpars,nh,its,ig)
      dimension par(kpars,nh,its),parj(nh,its)
     *         ,peion(4,nh,its)  !!! particle ionization  
      real i0,c(4)
      func1(e)=1.46e-15*(e**(-5.8)+5.e9*e**(-8.8))
!!    ! func2(e,e0,r)=e**(gam-0.854)*exp(-func1(e)*r**3-e/e0)
!!      data c/1.09,1.,0.61,1./
! NO particle ionization = 0
      data c/1.09,1.,0.61,0./

      a=1.277e-13*i0
      g=gam-0.854
      g1=1.+gam
      fg=gamma(g1)
      do 1 i=1,nh
        e1=emin
        r=1.09*par(1,i,ig)+par(2,i,ig)+0.61*par(3,i,ig)
	em=e0*(1.+8.16e-17*e0**0.44*r**1.323)
        em=0.373*em**0.1
        em=em*r**0.3
	em=em+g*e0

	psi=g+(5.76e-14*r)*r*r*(em**(-5.8)+1.1e10*em**(-8.8))
	
        sq=1./sqrt(psi)
        r=sqrt(psi/2.)
        er=1.+erf(r)
        fi=1.253/fg*sq*er
	fff=func1(em)
	fff=fff*r*r*r
        if(fff.gt.30.)fff=30.
        qer=em**(-0.854)*exp(-fff)
        eme0=em/e0
        if(eme0.ge.80.)eme0=80.
	  sum=(eme0)**g1
        sum=fi*qer*sum
      
        sumsum=exp(-eme0)
        sum=sum*sumsum
        if(sum.lt.1.e-30)sum=0.
    3   as=a*sum
c
        parj(i,ig)=0.
        do 4 j=1,4
          jj=j
          if(j.eq.3) jj=4
          if(j.eq.4) jj=3
          parj(i,ig)=parj(i,ig)+c(jj)*par(jj,i,ig)*as
! particle ionization calculates separately
!          par(j+12,i,ig)=par(j+12,i,ig)+c(jj)*par(jj,i,ig)*as
          peion(jj,i,ig)=c(jj)*par(jj,i,ig)*as
    4   continue
    1 continue
      return
      end

