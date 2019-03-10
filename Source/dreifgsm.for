      subroutine dreifGSM(nomsl,ntsl,nl,ht,tt,dolm,vdv,vdu,
     *       dtt,par,nr,par1,PAR2,park,ks,q,u,kpart,ddolgt,kdf,ldor,
     *       isp,pole,nv,rmaxt,vert,nh,ntr,rads,nin,mi,mast)
      dimension ntsl(nl),ht(nv),tt(nv),vdv(nv),vdu(nv),
     *       park(ks),q(nv,nl),u(nl),kdf(20),pole(ldor/4),
     *       vert(kpart,nv),rads(nh),par1(nr),q1(140),v1(2),u1(2),
     *       par(nr),PAR2(NR),mast(40)
      data re/6371.02e5/,pi/3.14159265/,dp/6.2831853/
      nt=ntsl(nomsl)
      rg=180./pi
c     print 900,mast(1)
  900 format(' dreifi:',i5)
      k1=nt/2
      k2=(nt+1)/2
      k3=k1+1
      do 1 k=1,nt
        h=ht(k)
        t=tt( k)
        r=re+h
        or=1./r
        ro=re*or
        ros=ro*ro
        ct=cos(t)
        st=sin(t)
        rst=or/st
        qq=ros*ct
        if(k.ne.1.and.k.ne.nt)go to 2
c  Расчет без учета переноса плазмы за счет зонального дрейфа
c        vv=dolm
c
c
        vv=dolm-vdv(k)*rst*dtt*rg
c

        if(vv.lt.0.)vv=vv+360.
        if(vv.ge.360.)vv=vv-360.
    2   continue
        if(k.ne.1)go to 3
c
cccc      u1(1)=u(nomsl)
cccc      u1(2)=u(nomsl)
c
c  Расчет без учета переноса плазмы за счет меридионального дрейфа
c           u1(1)=re*or*st*st
c           u1(2)=re*or*st*st
c
c
          u1(1)=re*or*(st-sqrt(1.+3.*ct*ct)*or*vdu(1)*dtt)*st
          u1(2)=re*or*(st-sqrt(1.+3.*ct*ct)*or*vdu(nt)*dtt)*st
c
c         ustate=re*or*st*st
c         print*,'ustate=',ustate,'u1(1)=',u1(1),'u1(2)=',u1(2)
    3   continue
        q1(k)=qq
        if(k.ne.nt)goto5
          v1(2)=vv
    5   continue
        if(k.ne.1)go to 4
          v1(1)=vv
    4   continue
        if(k.eq.1)uu=u1(1)
        if(k1.ne.k2)goto6
          if(k.eq.k3)uu=u1(2)
    6   continue
        if(k1.eq.k2) go to 7
          if(k.ne.k2) go to 7
            ht(k)=re*(1./uu-1.)
            tt(k)=pi*.5
            go to 1
    7   call njuan(qq,uu,h,t,rads(ntr))
        ht(k)=h
        tt(k)=t
    1 continue
      nfile=13
      if( nin.eq.0) nfile=6
c     print*,'dreif:intt: ntr,nr', ntr,nr
      call inttgsm (nomsl,vdu,vdv,q1,v1,u1,nfile,q,u,kpart,PAR,
     *       PAR2,ddolgt,rmaxt,ntsl,nl,kdf,ldor,isp,pole,ks,
     *       nv,vert,nh,ntr,park,rads,nr,PAR1,mi,ht,mast)
      return
      end

