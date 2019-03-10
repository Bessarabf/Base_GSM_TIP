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
