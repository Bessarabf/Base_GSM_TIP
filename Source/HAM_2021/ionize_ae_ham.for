!       partical ionization calculate separately 
!       massive peion 22/03/2019
!       add to interface ddolgs 10.07.2015  
      subroutine ionizE_AE_ham(par,parj,peion,ut,dolm,ddolgs,dtets,del,
     *                         nh,its,kpars,ps,E0,FAE)
      dimension ps(10),parj(nh,its),par(kpars,nh,its)
!                                  massive partical ionization for 4 species
     *         ,peion(4,nh,its)
      dimension E0(its,*),FAE(its,*)
      real ic,ia,ip,ifo

      data  emin/100./,emax/5.e4/,h/ 1./,ba/5./
      data   pi/3.141592/
! ps( ) defines in iacflo (iacflo_ae)
!  cusp presipitations
      ic=ps(1)   ! soft electron precipitation flux
      gamc=ps(2)
      e0c=ps(3)  ! characteristic energy
! auroral electrons 
      ia=ps(4)   ! auroral electrons
      gama=ps(5)
      e0a=ps(6)
! additional precipitations
      ip=ps(7)
      e0p=ps(8)
      ifo=ps(9)
! background precipitations
      e0f=ps(10)

      jg=dolm/ddolgs+1
c
      j=1
      phid=dolm*pi/180.
!  solar-magnetic coordinates
      call magsm(ut,del,phid,phism,j)

!!!! massive summary partical ionization
      parj=0.
     
!!!!      

      do 1 ig=1,its
        tet=dtets*(ig-1)

        ia=FAE(ig,jg)		   ! model Pakston
        e0a=E0(ig,jg)*1.e3
	
        qfsm=qf(ifo,phism,tet)
        if(tet.ge.90.-ba.and.tet.le.90.+ba) go to 2
c
          qcsm=qc(ic,phism,tet)
          qasm=ia

!       vertical part of ionization function
c         cusp precipitation (old var)
          call ionizv(qcsm,gamc,e0c,emin,emax,h,par,parj, 
     *                    kpars,nh,its,ig)

c         auroral precipitation separate (new var for EAGLE)
          if(e0a.ne.0.)call ionizv_ham(qasm,gama,e0a,emin,emax,h,par,
     *                             parj,peion,kpars,nh,its,ig)
           
    2   continue
    1 continue
      return
      end
