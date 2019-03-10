      subroutine insti(ntsl,nqi,par,pari,ni,park,kparn,rads,nh,ns,
     *ncs,dtets,kpars)
      dimension ntsl(nqi),par(kpars,nh,ncs), pari(ni),park(ns),
     *       rads(nh),msum(45),plm(8)
      msum (1)=0
      do 1 nomsl=2,nqi
        msum(nomsl)=msum(nomsl-1)+ntsl(nomsl-1)
    1 continue
      l=1
      do 6 nomsl=1,nqi
        nsum=msum(nomsl)*kparn
        nt=ntsl(nomsl)
        do 5 i=1,nt
          h=park(l)
          t=park(l+1)
          tr=abs(t-90.)
          if(tr.ge.0.01)go to 2
            call find(nh,h,rads,m)
            k=1
            ntet=90./dtets+1
            xi=h
            x1=rads(m)
            x2=rads(m+1)
            call int1(k,m,ntet,x1,x2,xi,par,plm,kparn,
     *                kpars,nh,ncs)
              go to 3
    2       continue
              m=i
              if(t.lt.90.)m=nt-i+1
              ntet=t/dtets
              xi=t
              x1=ntet*dtets
              x2=x1+dtets
              ntet=ntet+1
              k=2
              call int1(k,m,ntet,x1,x2,xi,par,plm,kparn,
     *                    kpars,nh,ncs)
    3       continue
            call vplm(plm(6),plm(7),t)
            ll=nsum+kparn*(i-1)+1
            do 4 j=1,kparn
              pari(ll)=plm(j)
              ll=ll+1
    4       continue
            l=l+2
    5   continue
    6 continue
      return
      end
