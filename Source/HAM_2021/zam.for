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

