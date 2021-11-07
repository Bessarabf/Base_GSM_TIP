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
