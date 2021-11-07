      subroutine inst(dolm,ntsl,nl,ntr,ddolgt,kdf,ldor,
     *       isp,par,pari,pole,nr,ni,park,ks,int,
     *       rads,nh,its,dtets,kpars,qom,qmax,iqo,mast)
      dimension ntsl(nl),par(kpars,nh,its), pari(ni),park(ks),
     *       rads(nh),pole(ldor/4),msum(45),plm(14),qom(nl)
     *       ,kdf(20),mast(40)
      logical readfl
  900 format(' ',10g12.4)
      msum (1)=0
      do 1 nomsl=2,nl
        msum(nomsl)=msum(nomsl-1)+ntsl(nomsl-1)
    1 continue
      l=1
      do 2 nomsl=1,nl
        nsum=msum(nomsl)*int
        nt=ntsl(nomsl)
        do 3 i=1,nt
         h=park(l)
         t=park(l+1)
         tr=abs(t-90.)
c  . . .  интерполяция ниже 520 км
         IF(h.le.rads(nh)) THEN
c  . . .    интерполяция на экваторе
            if(tr.lt.0.01) then
              call find(nh,h,rads,m)
              k=1
              ntet=90./dtets+1
              xi=h
              x1=rads(m)
              x2=rads(m+1)
              if(mast(32).eq.0) then
               call inter1(k,m,ntet,x1,x2,xi,par,plm,int,
     *                     kpars,nh,its)
              else
               call inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *                      kpars,nh,its)
              end if
            else
c  . . .     интерполяция ВНЕ экватора
              m=i+ntr-1
              if(t.lt.90.) m=nt-i+ntr
              ntet=t/dtets
              xi=t
              x1=ntet*dtets
              x2=x1+dtets
              ntet=ntet+1
              k=2
              if(mast(32).eq.0) then
	        
               call inter1(k,m,ntet,x1,x2,xi,par,plm,int,
     *                     kpars,nh,its)
              else
               call inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *                      kpars,nh,its)
              end if
            end if
          ELSE
c  . . .   интерполяция выше 520 км
           m=nh
           if(tr.lt.0.01) then
c  . . .    интерполяция на экваторе
             it=its/2+1
             do lm=1,3
               plm(lm)=par(lm,nh,it)
             end do
             plm(4)=par(6,nh,it)
             plm(5)=par(16,nh,it)
             plm(7)=2.e6
             plm(8)=par(7,nh,it)
             plm(14)=par(19,nh,it)
             if(mast(32).eq.0) then
               step=28.9*plm(8)**(-0.25)
               plm(6)=10.**step
c              plm(6)=(10.**step)*1.e-1
             else
               plm(6)=par(5,nh,it)
             end if
             plm(9)=par(10,nh,it)
             plm(10)=par(11,nh,it)
             plm(11)=par(12,nh,it)
             plm(12)=par(13,nh,it)+par(14,nh,it)+par(15,nh,it)
             call hplm (plm,rads(nh),h,int)
           else
c  . . .   интерполяция ВНЕ экватора
             ntet=t/dtets
             xi=t
             x1=ntet*dtets
             x2=x1+dtets
             ntet=ntet+1
             k=2
c  . . .    интерполяция, если Н по MSIS
             if(mast(32).eq.0) then
               call inter1(k,m,ntet,x1,x2,xi,par,plm,int,
     *                     kpars,nh,its)
             else
               call inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *                      kpars,nh,its)
             end if
             call hplm(plm,rads(nh),h,int)
           end if
         END IF
c    Сдвиг фазы меридионального ветра в Millstone Hill
cc         if(dolm.eq.0..or.dolm.eq.345.)then
cc           if(nomsl.eq.7.or.nomsl.eq.8)then
cc             plm(10)=plm(10)+6.e3
cc           end if
cc         end if
c
         call vplm (plm(9),plm(10),t)
         ll=nsum+int*(i-1)+1
         do 10 j=1,int
           pari(ll)=plm(j)
           ll=ll+1
   10    continue
         l=l+2
    3   continue
    2 continue
      j1=nl-1
      do13j=4,j1
        i1=0
        i2=j-1
        do12i=1,i2
          i1=i1+ntsl(i)
   12   continue
        i=iqo
        k=(i1+i-1)*int+1
        qo=pari(k+4)
        if(qom(j).lt.qo)qom(j)=qo
        if(qmax.lt.qo)qmax=qo
        n=ntsl(j)
        i=n-iqo+1
        k=(i1+i-1)*int+1
        qo=pari(k+4)
        if(qom(j).lt.qo)qom(j)=qo
        if(qmax.lt.qo)qmax=qo
   13 continue
      readfl=.false.
      nfile=8
      kpar=int
      md=1
      call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,
     * ntsl,nl,kdf,ldor,isp,md,pari,pole,ni,mast)
      return
      end
