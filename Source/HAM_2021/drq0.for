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

