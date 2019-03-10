      subroutine anals1(verno,god,day,ut0,utk,solen,dtt,dts,tau,
     *                 sole,solu,bmpz,bmpy,vsol,csol,fa0,fs,ap0,
     *                 pkp0,dst0,ae0,al0,au0,ntr,gamma,ddolgt,dtett,
     *                 nh,dh,rmin,b,c,nsu,nse,ddolgs,dtets,dlam,
     *                 nsill,dteta,ns,uprt,uprs,mas,mast,mass,rmaxt)
      dimension sole(nse),solu(nsu),mas(10),mast(40),mass(30)
      dimension solen(nse)
      integer uprt,uprs,verno,god,day
      integer mx(20)/1,1,1,100,20,15*0/,min(20)/20*0/
      integer mi(20)/20*0/,max(20)/3*20,1,20,1,2*20,1,20,1,20,
     *                    1,2,24,8,3,3*1/
c 902 format('  god  = ',i5)
c 903 format('  day  = ',i5)
c 905 format('   anals1   :  end   ')
c 906 format('  sole(',i3,') = ',e12.5)
c 907 format('  solu(',i2,')= ',e12.5)
c 908 format('  ut0    = ',f7.0)
c 909 format('  tk     = ',f7.0)
c 910 format('  utk    = ',f7.0)
c 911 format('  dts    = ',f7.0)
c 912 format('  dtt    = ',f7.0)
c 913 format('  fa0    = ',f7.0)
c 914 format('  fs     = ',f7.0)
c 915 format('  bmpy   = ',f7.0)
c 916 format('  bmpz   = ',f7.0)
c 917 format('  vsol   = ',e12.5)
c 918 format('  csol   = ',f7.0)
c 919 format('  ap0    = ',f7.0)
c 920 format('  pkp0   = ',f7.0)
c 921 format('  dst0   = ',f7.0)
c 922 format('  ae0    = ',f7.0)
c 923 format('  au0    = ',f7.0)
c 924 format('  al0    = ',f7.0)
c 925 format('  dlam   = ',f7.0)
c 926 format('  ddolgt = ',f7.0)
c 927 format('  ddolgs = ',f7.0)
c 928 format('  dlam   = ',f7.0)
c 929 format('  ntr    = ',i5)
c 930 format('  gamma  = ',e12.5)
c 931 format('  nh     = ',i5)
c 932 format('  rmin   = ',f5.0)
c 933 format('  dtett  = ',f5.0)
c 934 format('  dtets  = ',f5.0)
c 935 format('  dh     = ',f8.0)
c 936 format('  nsill  = ',i5)
c 937 format('  dteta  = ',f5.0)
c 938 format('  ns     = ',i5)
c 939 format('  uprt   = ',i3)
c 940 format('  uprs   = ',i3)
c 942 format('  l      = ',i3)
c 943 format('  mas(',i2,') = ',i3)
c 948 format('  uprt=',i2,' uprs=',i2,' if uprt=3 then  ',
c    *       'uprs=3')
c 949 format('  uprs=',i2,' uprt=',i2,' if uprs=3 then  ',
c    *       'uprt=3')
c 950 format('   anals1   :  begin ')
c 952 format('  tau=',f6.1)
c 953 format('  rmaxt  =',f8.0)
c 954 format('  rmint  =',f8.0)
c 955 format('  mast(',i2,') = ',i3)
c 956 format('  b',7x,'=',e12.5)
c 957 format('  c',7x,'=',e12.5)
c 958 format('  ddolgt.ne.dolgs !')
c 959 format('  solen(',i3,') = ',e12.5)
c 960 format('  for p/p ints:',
c    *       '   c-90.<dtets ! '/'  dtets=',e10.3,'   c=',e10.3)
c 961 format('  for p/p ints:',
c    *       '  180.-b<dtets ! '/'   dtets=',e10.3,'   b=',e10.3)
c 962 format('  incorrect  dtets=',e10.3)
c 963 format('  incorrect  ddolgt = ',e10.3)
  902 format('  anals1:  ******* error ******  god  = ',i5)
  903 format(' anals1: ******* error ******   day = ',i5)
  905 format('   anals1   :  end   ')
  906 format(' anals1: ******* error ******   sole(',i3,') = ',e12.5)
  907 format(' anals1: ******* error ******   solu(',i2,')= ',e12.5)
  908 format(' anals1: ******* error ******   ut0 = ',f7.0)
  909 format(' anals1: ******* error ******   tk = ',f7.0)
  910 format(' anals1: ******* error ******   utk = ',f7.0)
  911 format(' anals1: ******* error ******   dts = ',f7.0)
  912 format(' anals1: ******* error ******   dtt = ',f7.0)
  913 format(' anals1: ******* error ******   fa0 = ',f7.0)
  914 format(' anals1: ******* error ******   fs = ',f7.0)
  915 format(' anals1: ******* error ******   bmpy = ',f7.0)
  916 format(' anals1: ******* error ******   bmpz = ',f7.0)
  917 format(' anals1: ******* error ******   vsol = ',e12.5)
  918 format(' anals1: ******* error ******   csol = ',f7.0)
  919 format(' anals1: ******* error ******   ap0 = ',f7.0)
  920 format(' anals1: ******* error ******   pkp0 = ',f7.0)
  921 format(' anals1: ******* error ******   dst0 = ',f7.0)
  922 format(' anals1: ******* error ******   ae0 = ',f7.0)
  923 format(' anals1: ******* error ******   au0 = ',f7.0)
  924 format(' anals1: ******* error ******   al0 = ',f7.0)
  925 format(' anals1: ******* error ******   dlam = ',f7.0)
  926 format(' anals1: ******* error ******   ddolgt = ',f7.0)
  927 format(' anals1: ******* error ******   ddolgs = ',f7.0)
  928 format(' anals1: ******* error ******   dlam = ',f7.0)
  929 format(' anals1: ******* error ******   ntr = ',i5)
  930 format(' anals1: ******* error ******   gamma = ',e12.5)
  931 format(' anals1: ******* error ******   nh = ',i5)
  932 format(' anals1: ******* error ******   rmin = ',f5.0)
  933 format(' anals1: ******* error ******   dtett = ',f5.0)
  934 format(' anals1: ******* error ******   dtets = ',f5.0)
  935 format(' anals1: ******* error ******   dh = ',f8.0)
  936 format(' anals1: ******* error ******   nsill = ',i5)
  937 format(' anals1: ******* error ******   dteta = ',f5.0)
  938 format(' anals1: ******* error ******   ns = ',i5)
  939 format(' anals1: ******* error ******   uprt = ',i3)
  940 format(' anals1: ******* error ******   uprs = ',i3)
  942 format(' anals1: ******* error ******   l = ',i3)
  943 format(' anals1: ******* error ******   mas(',i2,') = ',i3)
  948 format(' anals1: ******* error ******   uprt=',i2,
     *       ' uprs=',i2,/' if uprt=3 then uprs=3')
  949 format(' anals1: ******* error ******   uprs=',i2,
     *       ' uprt=',i2,/' if uprs=3 then ', 'uprt=3')
  950 format('   anals1   :  begin ')
  952 format(' anals1: ******* error ******   tau=',f6.1)
  953 format(' anals1: ******* error ******   rmaxt =',f8.0)
  954 format(' anals1: ******* error ******   rmint =',f8.0)
  955 format(' anals1: ******* error ******   mast(',i2,') = ',i3)
  956 format(' anals1: ******* error ******   b',7x,'=',e12.5)
  957 format(' anals1: ******* error ******   c',7x,'=',e12.5)
  958 format(' anals1: ******* error ******   ddolgt.ne.dolgs !')
  959 format(' anals1: ******* error ******   solen(',i3,') = ',e12.5)
  960 format(' anals1: ******* error ******   for p/p ints:',
     *       '   c-90.<dtets ! '/'  dtets=',e10.3,'   c=',e10.3)
  961 format(' anals1: ******* error ******   for p/p ints:',
     *       '  180.-b<dtets ! '/'   dtets=',e10.3,'   b=',e10.3)
  962 format(' anals1: ******* error ******   dtets=',e10.3)
  963 format(' anals1: ******* error ******   ddolgt = ',e10.3)
      re=6371.02
      r15=15.*re
      print 950
      if(god.gt.1000.and.god.lt.2000) go to 1
      print 902,god
      verno=1
    1 continue
      if(day.gt.0.and.day.le.366)go to 2
      print 903,day
      verno=1
    2 continue
      m=mas(4)
      do 5 i=1,m
        if(sole(i).ge.0.)go to 5
      k=i
      print 906,k,sole(i)
      verno=1
    5 continue
      do 155 i=1,m
        if(solen(i).ge.0)go to 156
          print 959,i,solen(i)
          verno=1
  156   continue
  155 continue
      m=mas(5)
      do 8 i=1,m
        if(solu(i).ge.0.) go to 8
          k=i
          print 907,k,solu(i)
          verno=1
    8 continue
      m1=m+1
      m=m*2
      do 9 i=m1,m
        if(solu(i).ge.0.)go to 9
      k=i
      print 907,k,solu(i)
      verno=1
    9 continue
      if(ut0.ge.0..and.ut0.le.86400.)go to 20
      print 908,ut0
      verno=1
   20 continue
   21 if(utk.gt.0)go to 22
      print 910,utk
      verno=1
   22 if(dts.gt.0..and.dts.lt.7200.)go to 23
      print 911,dts
      verno=1
   23 if(dtt.gt.0..and.dtt.lt.7200.)go to 24
      print 912,dtt
      verno=1
   24 if(fa0.gt.0..and.fa0.lt.400.)go to 25
      print 913,fa0
      verno=1
   25 if(fs.gt.0..and.fs.lt.400.)go to 26
      print 914,fs
      verno=1
   26 continue
      if(bmpy.gt.-100..and.bmpy.lt.100.)go to 27
      print 915,bmpy
      verno=1
   27 if(bmpz.gt.-100..and.bmpz.lt.100.)go to 28
      print 916,bmpz
      verno=1
   28 if(vsol.gt.0.)go to 29
      print 917,vsol
      verno=1
   29 if(csol.lt.100..and.csol.gt.0.)go to 30
      print 918,csol
      verno=1
   30 if(ap0.gt.0..and.ap0.lt.1000.)go to 31
      print 919,ap0
      verno=1
   31 if(pkp0.gt.0..and.pkp0.lt.10.)go to 32
      print 920,pkp0
      verno=1
   32 if(dst0.lt.1000.)go to 33
      print 921,dst0
      verno=1
   33 if(ae0.lt.2 00..and.ae0.gt.0.)go to 34
      print 922,ae0
      verno=1
   34 if(au0.lt.2 00..and.au0.gt.0.)go to 35
      print 923,au0
      verno=1
   35 if(al0.lt.0..and.al0.gt.-1000.)go to 36
      print 924,al0
      verno=1
   36 continue
      if(dlam.gt.0..and.dlam.le.360.)go to 37
      print 925,dlam
      verno=1
   37 if(ddolgt.ge.15..and.ddolgt.le.90.)go to 38
      print 926,ddolgt
      verno=1
   38 if(ddolgs.ge.15..and.ddolgs.le.90.)go to 39
      print 927,ddolgs
      verno=1
   39 continue
      m=dlam
      mn=ddolgt
      mk=m/mn
      mn=m-mk*mn
      mc=0
      if(mn.eq.0)go to 40
      print 928,dlam
      mc=1
      verno=1
   40 continue
      m=dlam
      mn=ddolgs
      mk=m/mn
      mn=m-mk*mn
      if(mn.eq.0)go to 41
      if(mc.eq.1)go to 41
      print 928,dlam
      verno=1
   41 continue
      if(ntr.gt.0.and.ntr.le.119)go to 42
      print 929,ntr
      verno=1
   42 if(gamma.gt.0..and.gamma.lt.10.)go to 43
      print 930,gamma
      verno=1
   43 if(nh.gt.0.and.nh.le.30)go to 44
      print 931,nh
      verno=1
   44 if(dh  .gt.0..and.dh  .le.100.)go to 45
      print 935,dh
      verno=1
   45 if(rmin.gt.0..and.rmin.lt.200.)go to 46
      print 932,rmin
      verno=1
   46 if(dtett.ge.2..and.dtett.le.90.)go to 47
      print 933,dtett
      verno=1
   47 if(dtets.ge.5..and.dtets.le.45.)go to 48
      print 934,dtets
      verno=1
   48 if(nsill.gt.0.and.nsill.le.45)go to 49
      print 936,nsill
      verno=1
   49 if(dteta.gt.0..and.dteta.le.90.)go to 50
      print 937,dteta
      verno=1
   50 if(ns.le.10.and.ns.gt.0)go to 51
      print 938,ns
      verno=1
   51 if(uprt.eq.0.or.uprt.eq.1.or.uprt.eq.2.or.uprt.eq.3) go to 62
      print 939,uprt
      verno=1
      go to 63
   62 if(uprt.ne.3)go to 63
        if(uprs.eq.3)go to 64
      print 948,uprt,uprs
      verno=1
   64   continue
   63 if(uprs.eq.0.or.uprs.eq.1.or.uprs.eq.2.or.uprs.eq.3)go to 65
      print 940,uprs
      verno=1
      go to 66
   65 continue
      if(uprs.ne.3)go to 66
        if(uprt.eq.3)go to 167
      print 949,uprs,uprt
      verno=1
  167 continue
   66 continue
      mx(4)=nse
      mx(5)=nsu/2
      do 56 i=1,5
        if(mas(i).le.mx(i).and.mas(i).ge.mi(i))go to 57
          print 943,i,mas(i)
          verno=1
   57   continue
   56 continue
   60 if(tau.gt.0.)go to 70
      print 952,tau
        verno=1
   70 continue
      if(rmaxt.le.r15.and.rmaxt.gt.(re+200))go to 71
      print 953,rmaxt
      verno=1
   71 continue
      do 75 i=1,16
        if(mast(i).le.max(i).and.mast(i).ge.min(i))go to 76
        print 955,i,mast(i)
        verno=1
   76   continue
   75 continue
      b1=b
      if(b1.le.180..and.b1.gt.90.)go to 67
        print 956,b1
        verno=1
   67 continue
      c1=c
      if(c1.lt.180..and.c1.ge.90.)go to 68
        print 957,c1
        verno=1
   68 continue
      if(ddolgt.eq.ddolgs)go to 80
        print 958
        print 927,ddolgs
        print 926,ddolgt
        verno=1
   80 continue
      if(dtets.gt.(c-90.))go to 81
        print 960,dtets,c
        verno=1
   81 continue
      if(dtets.gt.(180.-b))go to 82
        print 961,dtets,b
        verno=1
   82 continue
      i=90./dtets
      i=90-i*dtets
      if(i.eq.0)go to 83
        print 962,dtets
        verno=1
   83 continue
      if(verno.eq.1)go to 61
        print 905
   61 continue
      return
      end
