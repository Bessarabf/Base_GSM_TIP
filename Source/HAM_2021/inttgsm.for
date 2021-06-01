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

