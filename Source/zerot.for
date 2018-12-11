      subroutine zerot(ddolgt,ntsl,nl,ut0,kdf,ldor,isp,kpart,
     *                 par,pole,nr,park,ks,mast)
       integer art
      logical readfl
      dimension kdf(20),par(nr),pole(ldor/4),msum(82)
      dimension ntsl(nl),park(ks),mast(40)
      real hmi(3)/300.,700.,500./
!      real cim(3)/5.e5,1.e1,1.e-3/
      real hi(3)/40.,400.,200./
      zbaz=125.
      tzbaz=300.
      tbesk=1000.
      s=0.014
      re=6371.02
      pi=3.1415926
      cr=pi/180.
      msum(1)=0
      do 10 i=2,nl
        msum(i)=msum(i-1)+ntsl(i-1)
   10 continue
cc    nfile=3
cc    kpar=2
cc    readfl=.true.
cc    dolm=0.
cc    md=1
cc    call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
cc   *       kdf,ldor,isp,
cc   *       md,park,pole,nr,mast)
      dolm=0.
    5 if(dolm.ge.360.)go to 1
        nomsl=1
   4    if(nomsl.gt.nl)go to 2
        l=kpart*msum(nomsl)
        l1=2*msum(nomsl)
        nt=ntsl(nomsl)
        do 3 i=1,nt
          m1=l1+2*(i-1)+1
          alt=park(m1) /1.e5
          m=l+kpart*(i-1)+1
          tet=park(m1+1)
          art=1
          call ggmraw(art,dol,tetg,dolm,tet)
          d=dol*cr+pi/2.
          pr=0.55+0.45*sin(d)
          do 6 k=1,3
            cnt=(alt-hmi(k))/hi(k)
            if(cnt.gt.70.)cnt=70.
c
ccc         par(m)=cim(k)*exp(1.-cnt-exp(-cnt))*pr
ccc         if(par(m).lt.1.e-25)par(m)=1.e-25
            par(m)=1.e-3
            m=m+1
    6     continue
          do 7 k=4,6
            par(m)=0.
            m=m+1
    7     continue
c
           sss=s*(alt-zbaz)
           if(sss.gt.70.)sss=70.
          tn=tbesk-(tbesk-tzbaz)*exp(-sss)
c         tn=tbesk-(tbesk-tzbaz)*exp(-s*(alt-zbaz))
          par(m)=tn
cc        par(m)=1000.
          par(m+1)=tn
cc        par(m+1)=1000.
    3   continue
        nomsl=nomsl+1
        go to 4
    2 continue
        kpar=kpart
        readfl=.false.
        md=1
        nfile=6
        call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
     *       kdf,ldor,isp,
     *       md,par,pole,nr,mast)
c
        dolm=dolm+ddolgt
        go to 5
    1 continue
      return
      end
