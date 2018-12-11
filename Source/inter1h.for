c . . . интерполяция для случая атом. Н по MSIS
      subroutine inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *        kpars,nh,its)
      dimension par(kpars,nh,its),pl(3),plm(int)
      if(k.eq.2) go to 1
        nt1=m
        nt2=m+1
        its1=ntet
        its2=its1
        go to 2
    1   continue
        its1=ntet
        its2=ntet+1
        nt1=m
        nt2=nt1
    2 continue
      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
      is=1
      j=1
      i=is
    9 if(i.gt.16)go to 3
c ******
      if(par(i,nt1,its1).le.0.or.par(i,nt2,its2).le.0.)print 909,
     *      i,nt1,its1,i,nt2,its2,par(i,nt1,its1),par(i,nt2,its2)
  909 format(' i,nt1,its1,i,nt2,its2,par(i,nt1,its1),par(i,nt2,its2)',
     *     ' - subr. inter1'/' ',6i4,2g12.3)
c*******
        pl(1)=alog(par(i,nt1,its1))
        pl(2)=alog(par(i,nt2,its2))
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(j)=exp(pl(3))
        if(i.eq.3) is=3
        if(i.eq.6) is=10
        i=i+is
        j=j+1
        go to 9
    3 continue
      plm(7)=2.e6
      is=3
      i=7
      j=8
   10 if(i.gt.12)go to 6
        pl(1)=par(i,nt1,its1)
        pl(2)=par(i,nt2,its2)
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(j)=pl(3)
        if(i.eq.10) is=1
        i=i+is
        j=j+1
        go to 10
    6 continue
c     step=28.9*plm(8)**(-0.25)
c     plm(6)=10.**step
c     plm(6)=(10.**step)*1.e-1
        pl(1)=alog(par(5,nt1,its1))
        pl(2)=alog(par(5,nt2,its2))
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(6)=exp(pl(3))
c
        argum=par(13,nt1,its1)+par(14,nt1,its1)+par(15,nt1,its1)
        pl(1)=alog(argum)
        argum=par(13,nt2,its2)+par(14,nt2,its2)+par(15,nt2,its2)
        pl(2)=alog(argum)
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(12)=exp(pl(3))
      return
      end

