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

