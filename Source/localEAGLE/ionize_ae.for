!       add to interface ddolgs 10.07.2015  
      subroutine ionize_AE(par,parj,ut,dolm,ddolgs,dtets,del,nh,its,
     *           kpars,ps,E0,FAE)
      dimension ps(10),parj(nh,its)
      dimension E0(its,*),FAE(its,*)
      real ic,ia,ip,ifo,par(kpars,nh,its)
      data  emin/100./,emax/5.e4/,h/ 1./,ba/5./
      data   pi/3.141592/
      ic=ps(1)   ! soft electron precipitation flux
      gamc=ps(2)
      e0c=ps(3)  ! characteristic energy
      ia=ps(4)   ! auroral electrons
      gama=ps(5)
      e0a=ps(6)
      ip=ps(7)
      e0p=ps(8)
      ifo=ps(9)
      e0f=ps(10)
!      jg=dolm/15.+1
      jg=dolm/ddolgs+1
c
      j=1
      phid=dolm*pi/180.
      call magsm(ut,del,phid,phism,j)
      do 1 ig=1,its
        tet=dtets*(ig-1)

        ia=FAE(ig,jg)		   ! model Pakston
        e0a=E0(ig,jg)*1.e3
	
        qfsm=qf(ifo,phism,tet)
        if(tet.ge.90.-ba.and.tet.le.90.+ba)go to 2
c
          qcsm=qc(ic,phism,tet)
          qasm=ia
c     !!! qp in south hemisphere is different!!!
          if(tet.le.90.) then
           qpsm=qp(ip,phism,tet)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

c         cusp precipitation
          call ionizv(qcsm,gamc,e0c,emin,emax,h,par,parj,kpars,
     *                nh,its,ig)

c         auroral precipitation
          if(e0a.ne.0.)call ionizv(qasm,gama,e0a,emin,emax,h,par,
     *                             parj,kpars,nh,its,ig)
    2   continue
c         ground precipitation
 !         call ionizv(qfsm,gamc,e0f,emin,emax,h,par,parj,kpars,
 !    *               nh,its,ig)
    1 continue
      return
      end
