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
        call spline(nto,x,y,b,c,d)
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

