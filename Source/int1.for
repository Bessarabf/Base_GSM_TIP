      subroutine int1(k,m,ntet,x1,x2,xi,par,plm,int,
     *        kpars,nh,its)
      dimension par(kpars,nh,its),pl(3),plm(int)
      if(k.eq.2) go to 1
        nt1=m
        nt2=m+1
        its1=ntet
        its2=its1
        go to 2
    1 continue
        its1=ntet
        its2=ntet+1
        nt1=m
        nt2=nt1
    2 continue
      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1.-w1
      is=1
      j=1
      i=is
    3 continue
      if(i.gt.6)goto 4
        if(par(i,nt1,its1).le.0.) print*,'int1 ',
     *                                    i,nt1,its1,par(i,nt1,its1)
        if(par(i,nt2,its2).le.0.) print*,'int1 ',
     *                                    i,nt2,its2,par(i,nt2,its2)
        pl(1)=alog(par(i,nt1,its1))
        pl(2)=alog(par(i,nt2,its2))
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(j)=exp(pl(3))
        if(i.eq.3) is=3
        i=i+is
        j=j+1
        go to 3
    4 continue
      is=1
      i=7
      j=5
    5 continue
      if(i.gt.12)go to 6
        pl(1)=par(i,nt1,its1)
        pl(2)=par(i,nt2,its2)
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(j)=pl(3)
        if(i.eq.7) is=3
        if(i.eq.10)is=1
        i=i+is
        j=j+1
        go to 5
    6 continue
      return
      end
