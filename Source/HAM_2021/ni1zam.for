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

