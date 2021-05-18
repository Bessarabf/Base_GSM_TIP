! subroutine dreifGSM, njuan, inttgsm, inttpl, inter4, inter6, inter8,
!      ni1raz, spline_bas, ni1zam, zamvv0, drq0, ni1net, zam, rastch
! function seval
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
!------------------------------------------------------------------------------
      subroutine njuan(qq,uu,oo,tt,h16)
      double precision x(4),q,u,c,s,r,p,v,a,b,apb,amb,o,t,rmv
      data re/6371.02e5/
      q=qq
      u=uu
      c=dble(1.)/dble(3.)
      s=u*u*.125/q
      r=s*s*.25
      p=dble(1.)/dble(27.)+r
      r=p+r
      v=dabs(s)*dsqrt(p)
      a=(r+v)**c
      rmv=r-v
      b=dabs(rmv)
      b=rmv/b*b**c
      v=a+b
      apb=v*.5
      amb=(a-b)*.5
      p=amb*dsqrt(dble(3.))
      r=c-apb
      a=datan2(p,r)*.5
      p=dsqrt(v+c)
      b=(amb*amb*3.+r*r)**.25*2.*dcos(a)
      apb=-p
      amb=-b
      x(1)=p+b
      x(2)=p+amb
      x(3)=apb+b
      x(4)=apb+amb
      h0=h16*1.e-5
      d0=0.d1
      do1i=1,4
        if(dabs(x(i)).gt.1.d0)goto1
          if(x(i).le.d0.and.s.gt.d0.or.x(i).ge.d0.and.s.lt.d0)goto1
            a=x(i)
            goto2
    1 continue
    2 continue
      s=dble(1.)-a*a
      p=dsqrt(s)
      t=datan2(p,a)
      o=u/s
      tt=t
      oo=re*(1./o-1.)
      return
      end
!------------------------------------------------------------------------------
!!!     VER 22.05.14 - par1 add to interface
!
      subroutine inttgsm(nomsl,vdu,vdv,q1,v1,u11,nfile,q,u,kpart,par,
     *           par1,ddolgt,rmaxt,ntsl,nl,kdf,ldor,isp,pole,ks,
     *        nv,vert,nh,ntr,park,rads,nr,par2,mi,ht1,mast)
!	include 'PARAMETR.INC'
      dimension vdu(nv),vdv(nv),q1(nv),v1(2),q(nv,nl),
     *          u(nl),park(ks),par(nr),par1(nr),par2(nr),ntsl(nl),
     *          kdf(20),u11(2),vert(kpart,nv),rads(nh),pole(ldor/4),
     *          hwx(2),ht1(nv),mast(40),msum(45)
      allocatable p1(:,:),p2(:,:),
     *            plm(:),plm1(:),plm2(:),plm3(:),
     *            plm4(:),plm5(:),plm6(:),plm7(:),
     *		    vn1(:),vn2(:),
     *            verh(:,:),x(:),y(:),b(:),c(:),d(:)
  
      data re/6371.02e5/

      allocate (p1(kpart,nv),p2(kpart,nv),
     *          plm(kpart),plm1(kpart),plm2(kpart),plm3(kpart),
     *          plm4(kpart),plm5(kpart),plm6(kpart),plm7(kpart),
     *		  vn1(kpart),vn2(kpart),verh(2,kpart),
     *          x(nv),y(nv),b(nv),c(nv),d(nv))
  
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iw=ntsl(nomsl)/2+1
      hhx=ht1(iw)
	msum(1)=0
      do 1 i=2,nl
        msum(i)=msum(i-1)+ntsl(i-1)
    1 continue
      nt=ntsl(nomsl)
      lr=0
      lin=1
      n1=1
      nn1=nt/2
      nn=nn1*2
      if(nn.ne.nt) go to 2
        n2=nn1
        j=2
        lin=0
        go to 5
    2 continue
        if(hhx.le.(rmaxt*1.e5-re)) go to 3
          n2=nn1+1
          j=2
          lr=1
          go to 4
    3   continue
          n2=nt
          j=1
    4   continue
    5 continue
      do 25 j1=1,j
        u1=u11(1)
        if(j1.eq.2)u1=u11(2)
c***
      hwx(1)=re*(1./u11(1)-1.)
      hwx(2)=re*(1./u11(2)-1.)
c***
        dolgx=v1(j1)
        ndp=dolgx/ddolgt
        dolg=ndp*ddolgt
        dolg1=dolg+ddolgt
        dolg2=dolg1
      if(abs(dolg1-360.).lt.0.1) dolg1=0.
        if(j1.eq.2) go to 6
c  Расчет без учета переноса плазмы за счет зонального дрейфа
c           vv=0.
c
c
          vv=vdv(1)
c
c  Расчет без учета переноса плазмы за счет меридионального дрейфа
c           vu=0.
          vu=vdu(1)
c
        go to 7
    6   continue
c  Расчет без учета переноса плазмы за счет зонального дрейфа
c           vv=0.
c
          vv=vdv(nt)
c
c  Расчет без учета переноса плазмы за счет меридионального дрейфа
c           vu=0.
c
          vu=vdu(nt)
c
    7   continue
        call inttpl(mi,dolg,dolg1,vv,par,par1,nr,
     *  nfile,kpart,ddolgt,ntsl,nl,kdf,POLE,ldor,isp,mast,dolga,dolgb)
        nob=msum(nomsl)*kpart
        if(vu.ne.0.) go to 17
          if(vv.ne.0.) go to 11
            if(j1.eq.2) go to 8
              l1=1
              go to 9
    8       continue
              l1=(n1-1)*kpart+1
    9       continue
            l2=n2*kpart
            do 10 i=l1,l2
              par2(nob+i)=par(nob+i)
   10       continue
            go to 15
   11     continue
            do 14 i=n1,n2
              m1=(i-1)*kpart
              do 12 m=1,kpart
                m2=m1+m+nob
                vn1(m)=par(m2)
                vn2(m)=par1(m2)
   12         continue
              call inter4(dolg,dolg2,dolgx,vn1,vn2,plm,kpart,MAST(24))
              do 13 m=1,kpart
                par2(nob+m1+m)=plm(m)
   13         continue
   14       continue
   15     continue
          go to 24
   17   continue
          if(u1.gt.u(1)) go to 18
            ni=0
            go to 19
   18     continue
            if(u1.le.u(nl)) go to 28
              ni=nl+1
              go to 19
   28       continue
            call find(nl,u1,u,ni)
   19     continue
          if(ni.ge.nl) go to 22
            ni1=ni+1
            nt=ntsl(ni1)
            nn=nt/2
            nn=nn*2
            if(nn.ne.nt) go to 20
            call ni1raz(ni,lr,n1,n2,ntsl,nl,vert,kpart,nv,ddolgt,
     *           nomsl,vv,p1,p2,par,par1,par2,nr,park,ks,vn1,vn2,
     *           plm,plm1,plm2,plm3,q1,u1,u,msum,dolg,dolg1,dolg2,
     *           dolgx,ht1,rmaxt,rads,nh,ntr,verh,x,y,b,c,d,MAST(24))
              go to 21
   20       continue
              call ni1zam(lin,j,lr,ni,n1,n2,ntsl,nl,vv,nomsl,u,u1,q1,
     *           plm1,plm2,plm3,plm4,plm,dolg,dolg1,dolg2,dolgx,ddolgt,
     *           msum,park,ks,par,par1,par2,nr,ht1,nv,kpart,p1,rads,nh,
     *           ntr,rmaxt,vn1,vn2,x,y,b,c,d,verh,plm5,plm6,plm7,
     *           MAST(24))
   21       continue
            go to 23
   22     continue
            call ni1net(ni,ntsl,nl,nomsl,n1,n2,rads,nh,ntr,
     *         u1,q1,dolg,dolg2,dolgx,vv,plm,plm1,plm2,plm3,
     *         msum,ht1,u,q,kpart,par,par1,par2,nr,nv,ddolgt,MAST(24))
   23     continue
   24   continue
        n1=n2
        if(lr.eq.0) n1=n2+1
        n2=ntsl(nomsl)
   25 continue
      if(j.ne.2) goto 27
        if(lr.ne.1) go to 29
          k=nn1*kpart
          do 26 m=1,kpart
            pl=(verh(1,m)+verh(2,m)) /2.
            par2(nob+k+m)=pl
   26     continue
   29   continue
   27 continue
  927 format(' ',8g12.3)
      
      deallocate (p1,p2,plm,plm1,plm2,plm3,
     *            plm4,plm5,plm6,plm7,
     *            vn1,vn2,verh,x,y,b,c,d)
      return
      end
!------------------------------------------------------------------------------
!      ver. 20.04.14 
	subroutine inttpl(mi,dolg,dolg1,vv,par,par1,nr,nfile,kpart,
     *           ddolgt,ntsl,nl,kdf,POLE,ldor,isp,mast,dolga,dolgb)
      dimension par(nr),par1(nr),ntsl(nl),kdf(20)
      dimension mast(40),pole(ldor/4)
      logical readfl
      readfl=.true.
        k1=3
        k2=3
        md=1
        if(mi.ne.1) go to 1
          call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *    kdf,ldor,isp,md,par,pole,nr,mast)
          dolga=dolg
          if(vv.eq.0.) go to 2
            call wwt(readfl,nfile,kpart,dolg1,ddolgt,1,ntsl,nl,
     *      kdf,ldor,isp,md,par1,pole,nr,mast)
            dolgb=dolg1
            go to 3
    2     continue
            do 4 i=1,nr
              par1(i)=par(i)
    4       continue
            dolgb=dolg
    3     continue
          mi=2
          go to 22
    1   continue
          if(vv.eq.0.) go to 8
            if(dolg1.eq.dolgb) go to 7
              k2=0
              if(dolg1.eq.dolga) k2=1
              dolgb=dolg1
    7       continue
            go to 9
    8     continue
            k2=2
    9     continue
          IF(dolg.NE.dolga) THEN
            k1=0
            if(dolg.eq.dolgb) k1=1
            dolga=dolg
          END IF
          if(vv.eq.0.) dolgb=dolg
          IF(k1.EQ.0) THEN
            if(k2.EQ.0) THEN
              call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par,pole,nr,mast)
              call wwt(readfl,nfile,kpart,dolg1,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par1,pole,nr,mast)
            end if
            IF(k2.EQ.1) then
              !!!!!!!!!!!!!!!!!!!do 13 i=1,nr
                par1=par
              !!!!!!!!!!!!!!!!!!!13         continue
              call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par,pole,nr,mast)
            END IF
            if(k2.EQ.2) THEN
              call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par,pole,nr,mast)
              !!!!!!!!!!!!!!!! do 14 i=1,nr
                par1= par !!!!!!!!!!par1(i)=par(i)
              !!!!!!!!!!!!!!!!!!14         continue
            END IF
           END IF
          if(k1.EQ.1) then
            if(k2.EQ.0) then
            !!!!!!!!!!!!!!!!!!!!!!          do 17 i=1,nr
              par=par1                     !!!par(i)=par1(i)
            !!!!!!!!!!!!!!!!!!!!!!!!!! 17         continue
              call wwt(readfl,nfile,kpart,dolg1,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par1,pole,nr,mast)
            end if  
            if(k2.eq.1) print 900
  900       format(' inttpl:    incorrect k2 ')
            if(k2.EQ.2) then
            !!!!!!!!!!!! do i=1,nr
                par=par1 !!!!!!!!!par(i)=par1(i)
            !!!!!!!!!!!!!19         continue
            end if
          end if
   22   continue
cc      mi=2
        return
        end
!------------------------------------------------------------------------------
      subroutine inter4(x1,x2,xi,vn1,vn2,plm,kpart,log)
      dimension plm(kpart),vn1(kpart),vn2(kpart)

      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
c
      if(log.eq.0)go to 11
        do 1 i=1,3
          p1=alog(vn1(i))
          p2=alog(vn2(i))
          plm(i)=exp(w1*p2+w2*p1)
    1   continue
        do 2 i=4,8
          plm(i)=w1*vn2(i)+w2*vn1(i)
    2   continue
        go to 12
   11 continue
        do 13 i=1,8
          plm(i)=w1*vn2(i)+w2*vn1(i)
   13   continue
   12 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine inter6(x1,x2,xi,vn1,vn2,kpart,log,nv,l,pm)
      dimension pm(kpart,nv),vn1(kpart),vn2(kpart)
      data al/0.54e-10/
      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
      do 1 i=1,3
        if(vn1(i).lt.al)vn1(i)=al
        if(vn2(i).lt.al)vn2(i)=al
        p1=vn1(i)
        p2=vn2(i)
        if(log.eq.0)go to 7
          p1=alog(vn1(i))
          p2=alog(vn2(i))
    7   continue
c       if(l.eq.5)print 77,l,nv,log,i,p1,p2
  77  format (' ',6g12.4)
        pm(i,l)=w1*p2+w2*p1
    1 continue
      do 2 i=4,8
        pm(i,l)=w1*vn2(i)+w2*vn1(i)
    2 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine inter8(plm1,plm2,plm3,kpart,uo,du,
     *                  dolo,ddolgt,ux,dolx,log,nv,l,pm,plm8)
	  
      dimension plm1(kpart),plm2(kpart),plm3(kpart),del(3)
      dimension pm(kpart,nv),plm8(kpart)
      data k/1/,al/0.54e-10/
      
      IF(ABS(DU).LT.(1.E-6))THEN
	      PRINT *,'INTER8 ',UO,DU,L,DOLO,DDOLGT,DOLX
!		 PAUSE
		 DU=1.E-6
		 dux=0.
      else
        dux=(ux-uo)/du
      END IF



      ddolx=(dolx-dolo)/ddolgt
      do 1 i=1,KPART !    15.04.15 ! I=1,8
        p1=plm1(i)
        p2=plm2(i)
        p3=plm3(i)
        p8=plm8(i)
        if(i.gt.3)go to 2
          if(p1.lt.al)p1=al
          if(p2.lt.al)p2=al
          if(p3.lt.al)p3=al
          if(p8.lt.al)p8=al
          if(log.eq.0)go to 3
            p1=alog(p1)
            p2=alog(p2)
            p3=alog(p3)
            p8=alog(p8)
   3      continue
   2    continue
        p4=p1+(p2-p1)*dux
        p5=p3+(p8-p3)*dux
        pm(i,l)=p4+(p5-p4)*ddolx

    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine ni1raz(ni,lr,n1,n2,ntsl,nl,vert,kpart,nv,ddolgt,
     *           nomsl,vv,p1,p2,par,par1,par2,nr,park,ks,vn1,vn2,
     *           plm,plm1,plm2,plm3,q1,u1,u,msum,dolg,dolg1,
     *           dolg2,dolgx,ht1,rmaxt,rads,nh,ntr,verh,
     *           x,y,b,c,d,log)
      dimension  ntsl(nl),vert(kpart,nv),p1(kpart,nv),
     *           p2(kpart,nv),par(nr),par1(nr),par2(nr),
     *           park(ks),vn1(kpart),vn2(kpart),plm(kpart),
     *           plm1(kpart),plm2(kpart),plm3(kpart),msum(45),
     *           u(nl),ht1(nv),y(nv),x(nv),q1(nv),b(nv),
     *           c(nv),d(nv),rads(nh),verh(2,kpart),plm8(8)
      double precision rr,tt,u1d,xd
      re=6371.02e5
      nn=ni
  900 format(' ni1raz')
5000  format(8(1pe10.2))
5001  format(' ')
c     print 900
      if(ni.eq.0)nn=1
      ii=ntsl(nn)
      nob1=msum(nn)
      nob2=msum(nn+1)
      if(n1.eq.1) go to 1
        i1=ii/2+1
        i2=ii
        go to 2
    1 continue
        i1=1
        i2=ii/2
    2 continue
      if(ni.ne.0) go to 30
        du=u(1)
        ul=0.
        go to 31
   30 continue
        du=u(ni+1)-u(ni)
        ul=u(ni)
   31 continue
      if(vv.eq.0.) go to 7
        do 10 in=i1,i2
          m1=(in-1)*kpart
          do 3 m=1,kpart
            m2=nob1*kpart+m1+m
            m3=nob2*kpart+m1+m
            if(ni.eq.0) go to 4
              plm1(m)=par(m2)
              plm2(m)=par(m3)
              go to 5
    4       continue
              plm1(m)=vert(m,in)
              plm2(m)=par(m2)
    5       continue
            plm3(m)=par1(m2)
            plm8(m)=par1(m3)
            if(ni.eq.0) plm3(m)=vert(m,in)
            if(ni.eq.0) plm8(m)=vert(m,in)
    3     continue
          call inter8(plm1,plm2,plm3,kpart,ul,du,dolg,
     *    ddolgt,u1,dolgx,log,nv,in,p1,plm8)
c  *** log
c          do 91 m=1,3
c            if(log.eq.1)p1(m,in)=alog(p1(m,in))
c  91      continue
   10   continue
        go to 6
    7 continue
        do 8 in=i1,i2
          m1=(in-1)*kpart
          do 24 m=1,kpart
            m2=nob1*kpart+m1+m
            m3 =nob2*kpart+m1+m
            if(ni.eq.0) go to 22
              plm1(m)=par(m2)
              plm2(m)=par(m3)
              go to 23
   22       continue
              plm1(m)=vert(m,in)
              plm2(m)=par(m2)
   23       continue
   24     continue
          call inter6(ul,u(ni+1),u1,plm1,plm2,kpart,log,nv,in,p1)
    8   continue
    6 continue
      nto=i2-i1+1
      ll=(i1-1)*2+nob1*2
      i=1
      j=1
      do 12 in=i1,i2
        rr=re/(re+park(ll+j))
        u1d=u1
        tt=dasin(dsqrt(u1d/rr))
        xd=rr*rr*dcos(tt)
        x(i)=xd
        if(n1.eq.1) x(i)=-x(i)
        i=i+1
        j=j+2
   12 continue
      nob1=msum(nomsl)*kpart
      do 11 m=1,kpart
        i=1
        do 13 in=i1,i2
          y(i)=p1(m,in)
          i=i+1
   13   continue
        call spline_bas(nto,x,y,b,c,d)
        do 14 n=n1,n2
          nn=nob1+(n-1)*kpart+m
          if(ht1(n).le.rads(ntr)) go to 15
            if(ht1(n).ge.(rmaxt*1.e5-re)) go to 16
              p=seval(nto,q1(n),x,y,b,c,d)
              go to 17
   16       continue
              ll=i2
              if(i1.ne.1) ll=i1
              p=p1(m,ll)
   17       continue
            go to 18
   15     continue
            ll=i1
            if(i1.ne.1) ll=i2
            p=p1(m,ll)
   18     continue
          if(m.gt.3) go to 271
            if(log.eq.0)go to 27
               par2(nn)=exp(p)
               go to 26
   27       continue
               par2(nn)=p
   26       continue
            go to 272
  271    continue
            par2(nn)=p
  272    continue
   14   continue
        if(lr.ne.1) go to 19
          ll=nob1+m+ntsl(nomsl)/2*kpart
          if(n1.ne.1) go to 20
            verh(1,m)=par2(ll)
            go to 21
   20     continue
            verh(2,m)=par2(ll)
   21     continue
   19   continue
   11 continue
c     n11=nob1+(n1-1)*kpart+m
c     n21=nob1+(n2-1)*kpart+m
c     print 5000,(par2(n),n=n11,n21)
c     print 5001
      return
      end
!------------------------------------------------------------------------------
      subroutine spline_bas(n,x,y,b,c,d)
      real x(n),y(n),b(n),c(n),d(n)
      data k/1/
      nm1=n-1
      if(n.lt.2) go to 999
      if(n.lt.3) go to 50
      d(1)=x(2)-x(1)
      c(2)=(y(2)-y(1))/d(1)
      do 10 i=2,nm1
        d(i)=x(i+1)-x(i)
        b(i)=2.*(d(i-1)+d(i))
        c(i+1)=(y(i+1)-y(i))/d(i)
        c(i)=c(i+1)-c(i)
   10 continue
      b(1)=-d(1)
      b(n)=-d(n-1)
      c(1)=0.
      c(n)=0.
      if(n.eq.3) go to 15
      c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1)=c(1)*d(1)**2/(x(4)-x(1))
      c(n)=-c(n)*d(n-1)**2/(x(n)-x(n-3))
   15 continue
      do 20 i=2,n
        t=d(i-1)/b(i-1)
        b(i)=b(i)-t*d(i-1)
        c(i)=c(i)-t*c(i-1)
   20 continue
      c(n)=c(n)/b(n)
      do 30 ib=1,nm1
        i=n-ib
        c(i)=(c(i)-d(i)*c(i+1))/b(i)
   30 continue
      b(n)=(y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.*c(n))
      do 40 i=1,nm1
        b(i)=(y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.*c(i))
        d(i)=(c(i+1)-c(i))/d(i)
        c(i)=3.*c(i)
   40 continue
      c(n)=3.*c(n)
      d(n)=d(n-1)
      go to 999
   50 continue
      b(1)=(y(2)-y(1))/(x(2)-x(1))
      c(1)=0.
      d(1)=0.
      b(2)=b(1)
      c(2)=0.
      d(2)=0.
  999 continue
      return
      end
!------------------------------------------------------------------------------
      real function seval(n,u,x,y,b,c,d)
      real x(n),y(n),b(n),c(n),d(n)
      data i/1/,k/1/
      if(i.ge.n)i=1
      if(u.lt.x(i)) go to 10
      if(u.le.x(i+1)) go to 30
   10 i=1
      j=n+1
   20 k=(i+j)/2
      if( u.lt.x(k) ) j=k
      if(u.ge.x(k)) i=k
      if(j.gt.(i+1)) go to 20
   30 dx=u-x(i)
      seval=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      k=k+1
      return
      end
!------------------------------------------------------------------------------
      subroutine ni1zam(lin,jc,lr,ni,n1,n2,ntsl,nl,vv,nomsl,u,u1,q1,
     *        plm1,plm2,plm3,plm4,plm,dolg,dolg1,dolg2,dolgx,ddolgt,
     *        msum,park,ks,par,par1,par2,nr,ht1,nv,kpart,p1,rads,nh,
     *        ntr,rmaxt,vn1,vn2,hi,y,b,c,d,verh,plm5,plm6,plm7,log)

      dimension ntsl(nl),plm1(kpart),plm2(kpart),plm3(kpart),
     *        plm(kpart),msum(nl),park(ks),par(nr),par1(nr),par2(nr),
     *        u(nl),ht1(nv),hi(nv),p1(kpart,nv),rads(nh),plm5(kpart),
     *        y(nv),b(nv),c(nv),d(nv),plm4(kpart),verh(2,kpart),
     *        vn1(kpart),vn2(kpart),q1(nv),plm6(kpart),plm7(kpart)
	allocatable qi(:)
 
      double precision rr,u1d,adwt,tt
      data alf/0.15/

	allocate (qi(nv))

      re=6371.02e5
      hhh=rmaxt*1.e5-re
      ni1=ni+1
      nob1=msum(ni)
      nob2=msum(ni1)
      nip=ntsl(ni)
      nip1=ntsl(ni1)
      if(n1.ne.1) go to 1
        i1=1
        if(jc.ne.2) go to 2
          i2=nip/2
          if(nip.ne.(i2*2)) i2=i2+1
          go to 3
    2   continue
          i2=nip
    3   continue
        go to 4
    1 continue
        i1=nip/2+1
        i2=nip
    4 continue
      ui=u(ni)
      ui1=u(ni1)
      du=ui1-ui
      mm=(msum(ni1)+ntsl(ni1)/2)*2+1
      hw1=park(mm)
      hwx=re*(1./u1-1.)
      k=nob1*2
      if(hwx.le.hhh)go to 5
        nus=nip/2
        nusd=nus
        go to 6
    5 continue
        nus=0
c*
        if(lin.eq.1)go to 8
          do 11 l=1,nip,2
            pp=park(k+l)
            if(pp.ge.hwx)go to 9
              nus=nus+1
    9       continue
   11     continue
          go to 12
    8   continue
c*
          do 40 l=1,nip,2
            pp=park(k+l)
c           if(pp.gt.hw1) go to 41
            if(pp.ge.hwx) go to 41
              nus=nus+1
   41       continue
   40     continue
   12   continue
c       if(lin.eq.1.and.hwx.le.(rmaxt*1.e5-re)) go to 7
c         nus=nus*2+1
c         nusd=nus/2+1
c         goto6
c   7   continue
        nus=nus*2+1
        nusd=nus/2+1
    6 continue
      l=1
      j=1
      if(i1.eq.1) go to 48
        j=i1*2-1
        nusd=1
   48 continue
      if(vv.ne.0.) go to 10
          call zamvv0(i1,i2,kpart,nob1,nob2,par,par1,nr,park,j,msum,
     *           nip,nip1,plm1,plm2,plm3,plm4,plm5,plm6,plm7,plm,ui,
     *           du,u1,u,dolg,dolg2,ddolgt,dolgx,l,nusd,jc,hwx,hw1,
     *        hi,ntsl,nomsl,nl,p1,vn1,vn2,nv,ks,ui1,vv,k,lin,log)
        go to 17
   10 continue
          call zam(i1,i2,kpart,nob1,nob2,par,par1,nr,
     *         park,j,msum,nip,nip1,plm1,plm2,plm3,plm4,plm5,plm6,
     *    plm7,plm,ui,du,u1,u,dolg,dolg2,ddolgt,dolgx,l,nusd,jc,hwx,
     *hw1,hi,ntsl,nomsl,nl,p1,vn1,vn2,nv,ks,k,re,lin,vv,log,rmaxt)
c********
      naa=dolgx/ddolgt+1
c     if(naa.eq.9.and.nomsl.eq.6)write(3,1000)(p1(1,iii),iii=1,l)
 1000 format(' ',7g12.4)
c****
   17 continue
      if(jc.eq.2) go to 101
        call drq0(u1,u,kpart,dolg,dolg2,dolgx,ddolgt,par,par1,nr,
     *  nomsl,ntsl,nl,vn1,vn2,plm1,plm2,plm3,vv,msum,nv,nusd,p1,log)
        indx=ntsl(nomsl)/2+1
        hi(nusd)=ht1(indx)
  101 continue
      nob1=msum(nomsl)*kpart
      u1d=u1
      do 27 i=1,l
        rr=re/(re+hi(i))
        adwt=dsqrt(u1d/rr)
        if(adwt.gt.0.1d1)adwt=0.1d1
        tt=dasin(adwt)
        qi(i)=rr*rr*dcos(tt)
        if(lin.eq.1.and.i.eq.nusd) qi(i)=0.
        if(jc.ne.2)go to 50
          if(n1.eq.1) qi(i)=-qi(i)
          go to 51
   50   continue
          lw=l/2+1
          if(i.lt.lw) qi(i)=-qi(i)
   51   continue
   27 continue
      do 32 m=1,kpart
        do 33 i=1,l
          y(i)=p1(m,i)
   33   continue
c****
c     if(naa.eq.9.and.nomsl.eq.6.and.m.eq.1)
c    * write(3,1000)l,(qi(iii),iii=1,l),(y(jjj),jjj=1,l)
c***
        call spline_bas(l,qi,y,b,c,d)
c****
c     if(naa.eq.9.and.nomsl.eq.6.and.m.eq.1)
c    * write(3,1000)l,n1,n2,(q1(iii),iii=n1,n2),(ht1(jjj),jjj=n1,n2),
c    *  rads(ntr),hi(nusd),nusd,(ht1(kkk),kkk=1,109)
        do 34 n=n1,n2
          nn=(n-1)*kpart
          nn=nob1+nn+m
          if(ht1(n).le.rads(ntr)) go to 35
            if(ht1(n).ge.hi(nusd)) go to 36
c****
c     if(naa.eq.9.and.nomsl.eq.6.and.m.eq.1.and.n.gt.54)
c    *  write(3,1000)l,n,q1(n),(qi(iii),iii=1,l),(y(jjj),jjj=1,l)
              p=seval(l,q1(n),qi,y,b,c,d)
c     if(naa.eq.9.and.nomsl.eq.6.and.m.eq.1.and.n.gt.54)
c    *  write(3,1000)l,n,q1(n),(qi(iii),iii=1,l),(y(jjj),jjj=1,l),p
              go to 37
   36       continue
              p=p1(m,nusd)
   37       continue
            go to 38
   35     continue
            ll=1
            if(i1.ne.1)ll=l
            p=p1(m,ll)
   38     continue
c   log |
          if(m.le.3.and.log.eq.1)p=exp(p)
          par2(nn)=p
   34   continue
c       kkk=msum(nomsl)*kpart+1
c       kkk1=kkk-1+ntsl(nomsl)*kpart
c       n3=n1+1
c       n4=n2-1
c       nn=(n3-1)*kpart
c       nn=nob1+nn+m
c       p=(1.-2.*alf)*par2(nn)+alf*(par2(nn+8)+par2(nn-8))
c       n3=n1+2
c       do 52 n=n3,n4
c          nn=(n-1)*kpart
c          nn=nob1+nn+m
c          pp=(1.-2.*alf)*par2(nn)+
c    *     alf*(par2(nn-8)+par2(nn+8))
c          par2(nn-8)=p
c          p=pp
c  52   continue
        if(lr.ne.1)go to 44
          ll=nob1+m+ntsl(nomsl)/2*kpart
          if(n1.ne.1)go to 45
            verh(1,m)=par2(ll)
            go to 46
   45     continue
            verh(2,m)=par2(ll)
   46     continue
   44   continue
   32 continue
   	deallocate (qi)
      return
      end
!------------------------------------------------------------------------------
      subroutine zamvv0(i1,i2,kpart,nob1,nob2,par,par1,nr,
     *    park,j,msum,nip,nip1,plm1,plm2,plm3,plm4,plm5,plm6,
     *    plm7,plm,ui,du,u1,u,dolg,dolg2,ddolgt,dolgx,l,nusd,jc,
     * hwx,hw1,hi,ntsl,nomsl,nl,p1,vn1,vn2,nv,ks,ui1,vv,k,lin,log)
      dimension ntsl(nl),plm1(kpart),plm2(kpart),plm3(kpart),plm(kpart),
     *   hi(nv),p1(kpart,nv),u(nl),msum(nl),park(ks),par(nr),par1(nr),
     *   vn1(kpart),vn2(kpart),plm5(kpart),plm6(kpart),plm7(kpart),
     *   plm4(kpart)
      data na/1/
      re=6371.02e5
        do 7 in=i1,i2
          h=park(k+j)
          m1=(in-1)*kpart
          m2=nob1*kpart+m1
          do 15 m=1,kpart
            plm1(m)=par(m2+m)
   15     continue
          if(h.gt.hw1) go to 8
            if(in.le.(nip1/2+1)) go to 64
              m1=(nip1-nip+in-1)*kpart
   64       continue
            m3=nob2*kpart+m1
            do 9 m=1,kpart
              plm2(m)=par(m3+m)
    9       continue
            call inter6(ui,ui1,u1,plm1,plm2,kpart,log,nv,l,p1)
            hi(l)=h
            l=l+1
            go to 11
    8     continue
            if(l.ne.nusd) go to 13
              if(jc.eq.2)go to 13
              l=l+1
              go to 14
   13       continue
c        if( jc.ne.2)go to 12
c             if(h.ge.hwx)go to 12
              if(h.gt.hwx)go to 12
                us=re/(re+h)
                vv0=0.
                call drq0(us,u,kpart,dolg,dolg2,dolgx,ddolgt,par,par1,
     *               nr,nomsl,ntsl,nl,vn1,vn2,plm5,plm6,plm7,vv0,msum,
     *               nv,l,p1,log)
                do 65 ij=1,kpart
                  plm4(ij)=p1(ij,l)
                  if(ij.le.3.and.log.eq.1)plm4(ij)=exp(p1(ij,l))
   65           continue
                call inter6(ui,us,u1,plm1,plm4,kpart,log,nv,l,p1)
                hi(l)=h
                l=l+1
   12        continue
   14       continue
   11     continue
          j=j+2
    7   continue
        l=l-1
c     do 93 ji=1,l
c         if(ji.eq.nusd.and.lin.eq.1)go to 94
c           do 82 ij=1,3
c             p1(ij,ji)=alog(p1(ij,ji))
c  82       continue
c  94     continue
c  93 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine drq0(ux,u,kpart,dolg,dolg2,dolgx,ddolgt,par,par1,nr,
     * nomsl,ntsl,nl,vn1,vn2,plm1,plm2,plm3,vv,msum,nv,lu,p1,log)
      dimension ntsl(nl),kdf(20),plm1(kpart),plm2(kpart),plm3(kpart),
     *p1(kpart,nv),vn 1(kpart),vn2(kpart),u(nl),msum(nl)
     *,par(nr),par1(nr),plm8(8)
      if(ux.ge.u(nl) ) go to 10
        do 1 k=1,nl
          l=(ntsl(k)/2)*2
          if(ntsl(k).eq.l) go to 2
            ni=k
            go to 55
    2     continue
    1   continue
   55   continue
        if(ux.ge.u(ni)) go to 9
          l=ntsl(ni)/2
          ni1=ni+1
          nob1=msum(ni)*kpart
          nob2=msum(ni1)*kpart
          l1=ntsl(ni1)/2
          mm1=l*kpart+nob1
          mm2=l1*kpart+nob2
          do 20 m=1,kpart
            vn2(m)=par(mm1+m)
            if(m.le.3)vn1(m)=1.e-3
            if(m.gt.3.and.m.lt.7)vn1(m)=0.
   20     continue
          vn1(7)=1000.
          vn1(8)=3000.
          un1=u(ni)
          un=1./14.9
          call inter6(un,un1,ux,vn1,vn2,kpart,log,nv,lu,p1)
          do 88 ij=1,kpart
             plm1(ij)=p1(ij,lu)
   88     continue
          if(vv.ne.0.) go to 3
            go to 5
    3     continue
          do 30 m=1,kpart
            vn2(m)=par1(mm1+m)
            if(m.le.3)vn1(m)=1.e-3
            if(m.gt.3.and.m.lt.7)vn1(m)=0.
   30    continue
          vn1(7)=1000.
          vn1(8)=3000.
          un1=u(ni)
          un=1./14.9
            call inter6(un,un1,ux,vn1,vn2,kpart,log,nv,lu,p1)
            do 31 lp=1,kpart
              plm2(lp)=p1(lp,lu)
              if(lp.le.3.and.log.eq.1)plm2(lp)=exp(plm2(lp))
              if(lp.le.3.and.log.eq.1)plm1(lp)=exp(plm1(lp))
   31       continue
            call inter6(dolg,dolg2,dolgx,plm1,plm2,kpart,
     *                     log,nv,lu,p1)
    5     continue
          go to 8
    9   continue
          call find( nl,ux,u,ni)
          l=ntsl(ni)/2
          ni1=ni+1
          l1=ntsl(ni1)/2
          mm1=(l+msum(ni) )*kpart
          mm2=(l1+msum(ni1))*kpart
          do 21 m=1,kpart
            plm1(m)=par(mm1+m)
            plm2(m)=par(mm2+m)
            plm3(m)=par1(mm1+m)
            plm8(m)=par1(mm2+m)
   21     continue
          un=u(ni)
          un1=u(ni1)
          if(vv.eq.0.) go to 11
            du=un1-un
            call inter8(plm1,plm2,plm3,kpart,un,du,dolg,
     *      ddolgt,ux,dolgx,1,nv,lu,p1,plm8)
            go to 6
   11     continue
            call inter6(un,un1,ux,plm1,plm2,kpart,1,nv,lu,p1)
    6     continue
    8   continue
        go to 7
   10 continue
        l=ntsl(nl)/2
        mm1=(l+msum(nl))*kpart
        do 22 m=1,kpart
          vn1(m)=par(mm1+m)
          vn2(m)=par1(mm1+m)
  22    continue
        if(vv.eq.0.)  go to 13
          call inter6(dolg,dolg2,dolgx,vn1,vn2,kpart,log,nv,lu,p1)
          go to 14
   13   continue
          do 12 m=1,kpart
            if(log.eq.0)go to 112
              if(m.le.3)go to 25
                p1(m,lu)=vn1(m)
                go to 26
   25         continue
                p1(m,lu)=alog(vn1(m))
   26         continue
              go to 113
  112       continue
              p1(m,lu)=vn1(m)
  113       continue
   12     continue
   14   continue
    7 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine ni1net(ni,ntsl,nl,nomsl,n1,n2,rads,nh,
     *           ntr,u1,q1,dolg,dolg2,dolgx,vv,plm,plm1,plm2,
     *           plm3,msum,ht1,u,q,kpart,par,par1,par2,nr,nv,ddolgt,log)
      dimension ntsl(nl),rads(nh),q1(nv),plm(kpart),plm1(kpart),
     *          plm2(kpart),plm3(kpart),msum(nl),ht1(nv),u(nl),q(nv,nl)
     *         ,par(nr),par1(nr),par2(nr)
     
      allocatable qq(:),plm8(:),pmn(:,:)
	allocate (qq(nv),plm8(kpart),pmn(kpart,nv))
      nob1=msum(nl)*kpart
      if(nomsl.ne.nl) go to 6
        if(vv.eq.0.) go to 3
          do 1 in=n1,n2
            mm=(in-1)*kpart+nob1
            do 2 m=1,kpart
              plm1(m)=par(mm+m)
              plm2(m)=par1(mm+m)
    2       continue
            call inter4(dolg,dolg2,dolgx,plm1,plm2,plm,kpart,log)
            do 4 m=1,kpart
              par2(mm+m)=plm(m)
    4       continue
    1     continue
          go to 33
    3   continue
          do 20 in=n1,n2
            mm=(in-1)*kpart+nob1
            do 21 m=1,kpart
              par2(mm+m)=par(mm+m)
   21       continue
   20     continue
   33   continue
        go to 5
    6 continue
        rr=rads(ntr)
        nn=ntsl(nomsl)
        n11=ntsl(nl)
        do 12 i=1,n11
          qq(i)=q(i,nl)
   12   continue
        nob2=msum(nomsl)*kpart
        nver=ntsl(nomsl)/2+1
        if(vv.eq.0.) go to 22
          do 8 in=1,nn
            mm1=(in-1)*kpart+nob2
            if(ht1(in).gt.rr) go to 7
              nq=1
              if(q1(in).gt.0.) nq=ntsl(nl)
              mm=(nq-1)*kpart+nob1
              do 9 m=1,kpart
                plm1(m)=par(mm+m)
                plm2(m)=par1(mm+m)
    9         continue
           call inter4(dolg,dolg2,dolgx,plm1,plm2,plm,kpart,log)
              do 11 m=1,kpart
                par2(mm1+m)=plm(m)
   11         con tinue
              go to 17
    7       continue
              qx=q1(in)
              if(in.ne.nver) go to 15
                mm=(in-1)*kpart+nob1
                do 31 m=1,kpart
                  plm1(m)=par(mm+m)
                  plm2(m)=par1(mm+m)
   31           continue
             call inter4(dolg,dolg2,dolgx,plm1,plm2,plm,kpart,log)
                do 555 m=1,kpart
                  par2(mm1+m)=plm(m)
  555           continue
                go to 17
   15         continue
                call find(n11,qx,qq,nq)
                mm=(nq-1)*kpart+nob1
                mm2=nq*kpart+nob1
                do 13 m=1,kpart
                  plm1(m)=par(mm+m)
                  plm2(m)=par(mm2+m)
                  plm3(m)=par1(mm+m)
                  plm8(m)=par1(mm2+m)
   13           continue
                q3=qq(nq)
                q2=qq(nq+1)
                dq=q2-q3
                call inter8(plm1,plm2,plm3,kpart,q3,dq,
     *          dolg,ddolgt,qx,dolgx,log,nv,in,pmn,plm8)
   30         continue
              do 16 m=1,kpart
                par2(mm1+m)=pmn(m,in)
                if(m.le.3)par2(mm1+m)=exp(par2(mm1+m))
   16         continue
   17       continue
    8     continue
          go to 23
   22   continue
          do 24 in=1,nn
            mm1=(in-1)*kpart+nob2
            if(ht1(in).gt.rr) go to 25
              nq=1
              if(q1(in).gt.0.) nq=ntsl(nl)
              mm=(nq-1)*kpart+nob1
              do 26 m=1,kpart
                par2(mm1+m)=par(mm+m)
   26         continue
              go to 27
   25       continue
              qx=q1(in)
              if(in.ne.nver) go to 10
                mm=(in-1)*kpart+nob1
                do 32 m=1,kpart
                  plm( m)=par(mm+m)
   32           continue
                go to 14
   10         continue
                call find(n11,qx,qq,nq)
                mm=(nq-1)*kpart+nob1
                mm2=nq*kpart+nob1
                do 28 m=1,kpart
                  plm1(m)=par(mm+m)
                  plm2(m)=par(mm2+m)
   28           continue
                q3=qq(nq)
                q2=qq(nq+1)
                call inter4(q3,q2,qx,plm1,plm2,plm,kpart,log)
   14         continue
              do 29 m=1,kpart
                par2(mm1+m)=plm(m)
   29         continue
   27       continue
   24     continue
   23   continue
    5 continue
      deallocate (qq,plm8,pmn)
      return
      end
!------------------------------------------------------------------------------
      subroutine zam(i1,i2,kpart,nob1,nob2,par,par1,nr,park,j,msum,nip,
     *nip1,plm1,plm2,plm3,plm4,plm5,plm6,plm7,plm,ui,du,u1,u,dolg,
     *dolg2,ddolgt,dolgx,l,nusd,jc,hwx,hw1,hi,ntsl,nomsl,nl,p1,vn1,vn2,
     *nv,ks,k,re,lin,vv,log,rmaxt)
      dimension ntsl(nl),plm1(kpart),plm2(kpart),plm3(kpart),
     *   msum(nl),park(ks),plm7(kpart),plm4(kpart),hi(nv),u(nl),
     *       vn2(kpart),plm(kpart),vn1(kpart),plm8(8),
     *   p1(kpart,nv),plm5(kpart),plm6(kpart),par1(nr),par(nr)
      do 50 iii=1,kpart
        p1(iii,nusd)=99.
   50 continue
        do 18 in=i1,i2
          m1=(in-1)*kpart
          m2=nob1*kpart+m1
          do 20 m=1,kpart
            plm1(m)=par(m2+m)
            plm3(m) =par1(m2+m)
   20     continue
          h=park(k+j)
          if(h.gt.hw1) go to 19
            if(in.le.(nip1/2+1)) go to 47
              m1=(nip1-nip+in-1)*kpart
   47       continue
            m3=nob2*kpart+m1
            do 21 m=1,kpart
              plm2(m)=par(m3+m)
              plm8(m)=par1(m3+m)
   21       continue
            call inter8(plm1,plm2,plm3,kpart,ui,du,dolg,
     *      ddolgt,u1,dolgx,log,nv,l,p1,plm8)
            hi(l)=h
            l=l+1
            go to 22
   19     continue
          if(l.ne.nusd) go to 24
            if(lin.ne.0)go to 23
              if(hwx.lt.(rmaxt*1.e5-re))go to 29
                 call rastch(re,h,ui,u,kpart,dolg,dolg2,dolgx,
     *           ddolgt,par,par1,nr,nomsl,ntsl,nl,vn1,vn2,
     * plm5,plm6,plm7,msum,nv,l,p1,plm4,plm8,plm1,plm3,u1,hi,log)
                 go to 30
   29          continue
                 call drq0(u1,u,kpart,dolg,dolg2,dolgx,ddolgt,par,par1,
     *                    nr,nomsl,ntsl,nl,vn1,vn2,plm1,plm2,plm3,vv,
     *                    msum,nv,nusd,p1,log)
                 hi(l)=hwx
                 l=l+1
   30          continue
               go to 25
   23        continue
               l=l+1
   25        continue
             go to 26
   24      continue
             if(h.ge.hwx)go to 27
               call rastch(re,h,ui,u,kpart,dolg,dolg2,dolgx,ddolgt,
     *         par,par1,nr,nomsl,ntsl,nl,vn1,vn2,plm5,plm6,plm7,
     *           msum,nv,l,p1,plm4,plm8,plm1,plm3,u1,hi,log)
   27        continue
   26      continue
   22     continue
          j=j+2
   18   continue
   17 continue
      l=l-1
      return
      end
!------------------------------------------------------------------------------
      subroutine rastch(re,h,ui,u,kpart,dolg,dolg2,dolgx,ddolgt,
     *         par,par1,nr,nomsl,ntsl,nl,vn1,vn2,plm5,plm6,plm7,
     *           msum,nv,l,p1,plm4,plm8,plm1,plm3,u1,hi,log)
      dimension ntsl(nl),plm1(kpart),plm3(kpart),
     *   msum(nl),plm7(kpart),plm4(kpart),hi(nv),u(nl),
     *   vn2(kpart),vn1(kpart),plm8(kpart),
     *   p1(kpart,nv),plm5(kpart),plm6(kpart),par1(nr),par(nr)
                us=re/(re+h)
                du1=us-ui
                v0=0.
                call drq0(us,u,kpart,dolg,dolg2,dolgx,ddolgt,par,par1,
     *          nr,nomsl,ntsl,nl,vn1,vn2,plm5,plm6,plm7,v0,msum,
     *          nv,l,p1,log)
                do 60 ij=1,kpart
                  plm4(ij)=p1(ij,l)
                  if(ij.le.3.and.log.eq.1)plm4(ij)=exp(plm4(ij))
c
   60           continue
                call drq0(us,u,kpart,dolg2,dolg2,dolgx,ddolgt,par1,par1,
     *          nr,nomsl,ntsl,nl,vn1,vn2,plm5,plm6,plm7,v0,msum,
     *          nv,l,p1,log)
                do 61 ij=1,kpart
                  plm8(ij)=p1(ij,l)
                  if(ij.le.3.and.log.eq.1)plm8(ij)=exp(plm8(ij))
c
   61           continue
                call inter8(plm1,plm4,plm3,kpart,ui,du1,dolg,
     *          ddolgt,u1,dolgx,log,nv,l,p1,plm8)
                hi(l)=h
                l=l+1
      return
      end