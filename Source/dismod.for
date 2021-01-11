      function dismod(ano,tem,g,rh,solu,nsu,hi,key)
      dimension solu(nsu),sp(12)
      data re/6.371e8/,am1/53.12e-24/,bk/1.38e-16/,
     *     pi/3.14159265359d0/,
cc      cross section for F10.7= 70
cc   *     sp/0.12e-18,1.39e-18,13.2e-18,10.5e-18,2.39e-18,
cc   *        0.87e-18,0.28e-18,0.01e-18,0.,0.,0.,0./
c       cross section for F10.7=115
     *     sp/7.29E-19,4.30E-18,4.03E-19,4.69E-19,2.29E-18,9.40E-18,
     *        1.37E-17,1.01E-17,5.98E-18,2.55E-18,1.08E-18,3.93E-19/
c       cross section for F10.7=180
cc   *     sp/8.14e-19,4.27e-18,4.51e-19,4.69e-19,2.31e-18,9.48e-18,
cc   *        1.37e-17,1.01e-17,6.01e-18,2.58e-18,1.09e-18,3.92e-19/
! cross section and spectral interval from Ackerman et al. Planet Space Sci. 1970. v. 1970 (20 intervals)
! and reduction to 12 intervals. Cross-sections averaging with flux value
!
! 1210.0  1220.0  
! 1220.0  1250.0
! 1250.0  1270.0
! 1270.0  1310.0
! 1310.0  1350.0  
! 1350.0  1380.0
! 1380.0  1500.0 
! 1500.0  1550.0
! 1550.0  1630.0 
! 1630.0  1670.0
! 1670.0  1720.0
! 1720.0  1760.0
      
      dis mod=0.
      sum=0.
!!! O2 scale height
      h=bk*tem/(am1*g)
!!!
      ra=sqrt(rh*(rh+2.*re))
      alfa=atan(re/ra)
      him=pi-alfa
      if((hi.gt.him)) dis mod=0.
      
      f=cos(hi)
c      arad=80.0/180.*pi
      arad=85.0/180.*pi
      if(.not.(hi.gt.arad)) go to 2
        z=(rh+re)*(1.-sin(hi))/h
        if(z.gt. 90.) z= 90.
        er=erf(sqrt(z))
        if(er.eq.1..and.f.gt.0.) go to 2
        sec=sqrt(pi/2.*re/h)*exp(z)
        sec=sec*(1.-sign(1.,f)*er)
        go to 5
    2 continue
      sec=1./f
    5 continue
      tt=ano*h*sec
      lin=nsu/key
      l0=1
      if(key.eq.1) l0=l0+nsu/2
      l1=1
      do 6 l=l0,lin
        tau=sp(l1)*tt
        sum=0.
        if(tau.le.60. )
     *    sum=sp(l1)*solu(l)*exp(-tau)
        dismod=dismod+sum
        l1=l1+1
    6 continue
      return
      end
c
