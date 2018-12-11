      subroutine g52t41(ni,nim,nj,x,y,jww,psdv)
      dimension x(nim),y(nim),psdv(nim,nim),pr(50,50)
      jww=52041
      mi=ni+1
      mj=nj+1
      do 2 j=1,nj
        do 1 i=1,ni
          pr(i,j)=psdv(i,j)
    1   continue
        pr(mi,j)=pr(1,j)
    2 continue
      do 3 i=1,ni
        pr(i,mj)=pr(i,1)
    3 continue
      pr(mi,mj)=pr(2,2)
      do 5 j=2,nj-1
        tg=y(j)
        do 4 i=1,ni
          dg=x(i)
          psdv(i,j)=g52t51(dg,ni,nim,nj,pr,tg,x,y)
    4   continue
    5 continue
      tg=0.
      dg=30.
      p=g52t51(dg,ni,nim,nj,pr,tg,x,y)
      do 6 i=1,ni
        psdv(i,1)=p
    6 continue
      tg=180.
      p=g52t51(dg,ni,nim,nj,pr,tg,x,y)
      do 7 i=1,ni
        psdv(i,nj)=p
    7 continue
      return
      end
      function g52t51(dg,ni,nim,nj,pr,tg,x,y)
      dimension pr(nim,nim),x(nim),y(nim)
      ks=0
      call g52t61(ks,dg,tg,dm,tm)
      i=j52t62(nim,ni,dm,x)
      j=j52t62(nim,nj,tm,y)
      k=i+1
      l=j+1
      r=x(i)
      cx=(dm-r)/(x(k)-r)
      r=y(j)
      cy=(tm-r)/(y(l)-r)
      r=1.-cy
      s=1.-cx
      a=r*s
      b=r*cx
      c=cy*s
      d=cy*cx
      r=pr(i,j)*a+pr(k,j)*b
      s=pr(i,l)*c+pr(k,l)*d
      g52t51=r+s
      return
      end
c
c             dolg,dolm - geograf.,geomagn. dolg
c             tet,tetm - geograf.,geomagn. colatidute
c           art=0 :      g --> m
c
      subroutine g52t61(art,dolg,tet,dolm,tetm)
      integer art
      double precision zpi,faktor,cbg,ci,si,xlm,bm,cbm,sbm,
     *                 clm,slm,sbg,bg,slg,clg,xlg,ylg
      zpi=6.28318530718
      faktor=0.01745329252
      cbg=11.4*faktor
      ci=dcos(cbg)
      si=dsin(cbg)
      if(art.eq.0) go to 10
        xlm=dolm
        bm=90.-tetm
        cbm=dcos(bm*faktor)
        sbm=dsin(bm*faktor)
        clm=dcos(xlm*faktor)
        slm=dsin(xlm*faktor)
        sbg=sbm*ci-cbm*clm*si
        bg=dasin(sbg)
        cbg=dcos(bg)
        slg=(cbm*slm)/cbg
        clg=(sbm*si+cbm*clm*ci)/cbg
        if(clg.gt.1..or.(1.-clg).lt.1.e-10)clg=1.
        if(clg.lt.-1..or.(1.+clg).lt.1.e-10)clg=-1.
        xlg=dacos(clg)
        if(slg.lt.0.0) xlg=zpi-dacos(clg)
        bg=bg/faktor
        xlg=xlg/faktor
        xlg=xlg-69.8
        if(xlg.lt.0.0) xlg=xlg+360.0
        tet=90.-bg
        dolg=xlg
        go to 20
   10   bg=90.-tet
        xlg=dolg
        ylg=xlg+69.8
        cbg=dcos(bg*faktor)
        sbg=dsin(bg*faktor)
        clg=dcos(ylg*faktor)
        slg=dsin(ylg*faktor)
        sbm=sbg*ci+cbg*clg*si
        bm=dasin(sbm)
        cbm=dcos(bm)
        slm=(cbg*slg)/cbm
        clm=(-sbg*si+cbg*clg*ci)/cbm
        xlm=dacos(clm)
        if(slm.lt.0.0) xlm=zpi-dacos(clm)
        bm=bm/faktor
        xlm=xlm/faktor
        dolm=xlm
        if(abs(tet-180.).lt.1.e-3)dolm=0.
        tetm=90.-bm
   20 continue
      return
      end
      function j52t62(m,n,u,x)
      dimension x(m)
      i=n
      s=x(i)
      i=n-1
      if(u.ge.s) goto 2
      i=1
      s=x(i+1)
      if(u.lt.s) goto 2
      j=n+1
    1 continue
        k=(i+j)/2
        s=x(k)
        if(u.lt.s) j=k
        if(u.ge.s) i=k
      if(j.gt.i+1) goto 1
    2 continue
      j52t62=i
      return
      end
