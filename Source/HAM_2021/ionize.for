      subroutine ionize(par,parj,ut,dolm,dtets,del,nh,its,kpars,ps)
      dimension ps(10),parj(nh,its)
      real ic,ia,ip,ifo,par(kpars,nh,its)
      data  emin/100./,emax/5.e4/,h/ 1./,ba/5./
      data   pi/3.141592/
! ps( ) defines in iacflo (iacflo_ae)
!  cusp presipitations
      ic=ps(1)     ! electron precipitation flux 
      gamc=ps(2)   ! gamma  
      e0c=ps(3)    ! characteristic energy
! auroral electrons 
      ia=ps(4)     ! 
      gama=ps(5)
      e0a=ps(6)
! additional precipitations
      ip=ps(7)
      e0p=ps(8)
      ifo=ps(9)
! background precipitations
      e0f=ps(10)
c
      j=1
      phid=dolm*pi/180.
      call magsm(ut,del,phid,phism,j)
      do 1 ig=1,its
        tet=dtets*(ig-1)
          qfsm=qf(ifo,phism,tet)
        if(tet.ge.90.-ba.and.tet.le.90.+ba)go to 2
c
          qcsm=qc(ic,phism,tet)  ! latitude and longitude dependence
          qasm=qa(ia,phism,tet)  !
c     !!! qp in south hemisphere is different!!!
          if(tet.le.90.) then
c          qpsm=0.
           qpsm=qp(ip,phism,tet)
           goto 3
          end if
           qpsm=qp(ip,phism,tet)
          if(dolm.ge.180.) then
            dolms=360.-dolm
          else
            dolms=dolm
          end if
          qpsm=qpsm*exp(-(dolms-10.)/30.)**2
    3   continue
c
c         cusp precipitation
          call ionizv(qcsm,gamc,e0c,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)

c         auroral precipitation
          call ionizv(qasm,gama,e0a,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)
    2   continue
    1 continue
c 667 format('667 - ionize')
c     print 667
      return
      end
