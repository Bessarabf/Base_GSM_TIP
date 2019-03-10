      subroutine ionize(par,parj,ut,dolm,dtets,del,nh,its,kpars,ps)
      dimension ps(10),parj(nh,its)
      real ic,ia,ip,ifo,par(kpars,nh,its)
      data  emin/100./,emax/5.e4/,h/ 1./,ba/5./
      data   pi/3.141592/
      ic=ps(1)
      gamc=ps(2)
      e0c=ps(3)
      ia=ps(4)
      gama=ps(5)
      e0a=ps(6)
      ip=ps(7)
      e0p=ps(8)
      ifo=ps(9)
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
          qcsm=qc(ic,phism,tet)
          qasm=qa(ia,phism,tet)
c     !!! qp in south hemisphere is different!!!
          if(tet.le.90.) then
c          qpsm=0.
           qpsm=qp(ip,phism,tet)
c          if(ut.ge.0..and.ut.le.25200.)then
c            qpsm=qpsm*.25*(5.+3.*sin(pi/46800.*ut-pi/26.))
c            qpsm=qpsm*2.
c          else
c            if(ut.gt.25200..and.ut.le.50400.)then
c              qpsm=qpsm*.25*(5.+3.*sin(pi/25200.*ut-pi*.5))
c            else
c              if(ut.ge.64800..and.ut.le.86400.)then
c              if(ut.ge.75600..and.ut.le.86400.)then
c                qpsm=qpsm*.25*(5.+3.*sin(pi/46800.*ut-49.*pi/26.))
c                qpsm=qpsm*3.
c              else
c                qpsm=qpsm*.5
c              end if
c            end if
c          end if
           goto 3
          end if
c         if(dolm.ge.315.or.dolm.le.60.)
c    *    qpsm=qp(ip,phism,tet)
           qpsm=qp(ip,phism,tet)
          if(dolm.ge.180.) then
            dolms=360.-dolm
          else
            dolms=dolm
          end if
          qpsm=qpsm*exp(-(dolms-10.)/30.)**2
c         qpsm=0.
    3   continue
c
c         ground precipitation
c         call ionizv(qfsm,gamc,e0f,emin,emax,h,par,parj,kpars,
c    *               nh,its,ig)
c         cusp precipitation
          call ionizv(qcsm,gamc,e0c,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)
c         addition precipitation
c         if(tet.ge.90.) then
c           call ionizv(qpsm,gamc,e0p,emin,emax,h,par,parj,kpars,
c    *               nh,its,ig)
c           else
c           !!different Eo=e0c for northen hemisphere!!
c           call ionizv(qpsm,gamc,e0c,emin,emax,h,par,parj,kpars,
c    *               nh,its,ig)
c         end if
c         auroral precipitation
          call ionizv(qasm,gama,e0a,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)
    2   continue
c         ground precipitation
          call ionizv(qfsm,gamc,e0f,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)
    1 continue
c 667 format('667 - ionize')
c     print 667
      return
      end
