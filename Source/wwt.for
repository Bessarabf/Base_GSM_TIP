      subroutine wwt (readfl,nfile1,kpar,dolm,ddolgt,nomsl,
     *       ntsl,nl,kdf,ldor,isp,md,par,pole,nr,mast)
      dimension ntsl(nl),kdf(20),par(nr),pole(ldor/4),mast(40)
      logical readfl
      nfile=nfile1
      if(nfile.eq.3)print 101,nfile,readfl
  101 format(' wwt :   nfile=',i3,'  readfl=',l8)
      if(nfile.eq.14)go to 8
      go to 9
    8 continue
      if(mast(13).eq.0) nfile=13
c     if(mast(1).eq.1)nfile=13
    9 continue
      nob1=0
      nob2=0
      ns=0
      do 1 i=1,nl
        ns=ns+ntsl(i)
    1 continue
      if(nfile.eq.3)go to 3
        if(dolm.eq.0.)go to 2
          nn=dolm/ddolgt
          nob2=kpar*ns*nn
    2   continue
    3 continue
      if(md.ne.0)go to 6
        if(nomsl.eq.1)go to 5
          ns1=0
          k=nomsl-1
          do 4 i=1,k
            ns1=ns1+ntsl(i)
    4     continue
          nob1=kpar*ns1
    5   continue
        lpar=ntsl(nomsl)*kpar
        go to 7
    6 continue
        lpar=kpar*ns
    7 continue
      nob=nob1+nob2
      call inpout(readfl,nfile,kdf,ldor,
     *       isp,nob,lpar,par,pole,nr)
      return
      end

