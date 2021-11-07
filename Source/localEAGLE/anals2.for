      subroutine  anals2(kpart,kpars,ddolgt,ddolgs,
     *dtets,nh,ntsl,nl,dteta,dlam,nsill,
     *kdf,ldor,isp,par,pole,nr,verno,pari,nv,its,mast,mass)
      integer verno
      logical readfl
      dimension pari(kpart,nv),mast(40),mass(30)
      dimension ntsl(nl),kdf(20),par(kpars,nh,its),pole(ldor/4)
  899 format('   anals2   :  begin    ')
  901 format('  par(',i2,',',i2,',',i2,')=',e12.5)
  902 format('  par(',i2,',',i3,')=',e12.5,5x,' nomsl=',i2)
  903 format('  dolg=',e12.5)
  111 format('   anals2   :  incorrect   ,    STOP ! ')
  999 format('   anals2   :  end  ')
      print 899
      readfl=.true.
      verno=0
      skorn=1.e6
      skori=1.e20
      coni =1.e8
      tempn=1.e4
      tempe=1.e6
      tempi=1.e5
      qis=1.e5
      nfile=5
      idm=360./ddolgs
      idolg=dlam/ddolgs
      itet=dteta/dtets
      md=1
      iv1=0
      dolg=0.
      do 102 id=1,idm,idolg
        call wws(readfl,nfile,kpars,dolg,ddolgs,tet,dtets,nh,
     *  kdf,ldor,isp,md,par,pole,nr)
        do 101 it=1,its,itet
          do 100 l=1,nh
            do 10 k=1,kpars
              a=par(k,l,it)
              if(k.gt.6) go to 11
                if(a.lt.0.) verno=1
   11         continue
              if(k.ne.7) go to 12
                if(a.lt.10..or.a.ge.tempn) verno=1
c               if(a.lt.0..or.a.ge.tempn) verno=1
   12         continue
              if(k.ne.8) go to 13
                if(a.lt.10..or.a.ge.tempi) verno=1
c               if(a.lt.0..or.a.ge.tempi) verno=1
   13         continue
              if(k.ne.9) go to 14
                if(a.lt.10..or.a.ge.tempe) verno=1
c               if(a.lt.0..or.a.ge.tempe) verno=1
   14         continue
              if(k.lt.10.or.k.gt.12) go to 15
                if(abs(a).ge.skorn) verno=1
   15         continue
              if(k.gt.16.or.k.lt.13) go to 16
                if(a.lt.0..or.a.ge.qis) verno=1
   16         continue
              if(verno.eq.1) print 901,k,l,it,a
   10       continue
            if(verno.eq.0) go to 17
              verno=0
              iv1=1
              print 903,dolg
   17       continue
  100     continue
  101   continue
        dolg=dolg+dlam
  102 continue
      nfile=6
      md=0
      idm=360./ddolgt
      idolg=dlam/ddolgt
      dolg=0.
      do 103 id=1,idm,idolg
        do 104 nomsl=1,nl,nsill
        call wwt(readfl,nfile,kpart,dolg,ddolgt,
     *    nomsl,ntsl,nl,kdf,ldor,isp,
     *    md,pari,pole,nr,mast)
          nt=ntsl(nomsl)
          do 1 l=1,nt
            do 2 k=1,kpart
              a=pari(k,l)
              if(k.gt.3) go to 3
                if(a.lt.0..or.a.ge.coni) verno=1
    3         continue
              if(k.lt.4.or.k.gt.6) go to 4
                if(abs(a).ge.skori) verno=1
    4         continue
              if(k.ne.7) go to 5
                if(a.lt.0..or.a.ge.tempi) verno=1
    5         continue
              if(k .ne.8) go to 6
                if(a.lt.0..or.a.ge.tempe) verno=1
    6         continue
              if(verno.eq.1) print 902,k,l,a,nomsl
    2       continue
            if(verno.eq.0) go to 7
              print 903,dolg
    7       continue
    1     continue
  104   continue
        dolg=dolg+dlam
  103 continue
      if(iv1.eq.1)verno=1
      if(verno.eq.1) print 111
      print 999
      return
      end

