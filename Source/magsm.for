      subroutine magsm(ut,del,phid,phism,j)
      data g10/30103.6/,g11/-2016.5/,h11/5682.6/,pi2/6.283185/
  900 format(' ',10g12.4)
      sq=g11**2+h11**2
      sqq=sqrt(sq)
      sq=sqrt(sq+g10**2)
      st0=sqq/sq
      ct0=g10/sq
      sl0=-h11/sqq
      cl0=-g11/sqq
      hour=ut/3600.
      al=0.2618*(hour-12.)
      sal=sin(al)
      cal=cos(al)
      ssd=sin(del)
      csd=cos(del)
      x1=csd*cal
      y1=-csd*sal
      z1=(x1*cl0+y1*sl0)*ct0-ssd*st0
      y1=y1*cl0-x1*sl0
      xmut=12.-3.8197*atan2(y1,z1)
      fi=-0.2618*xmut
      if(j.lt.0) go to 1
      phism=phid-fi-pi2/2.
      if(phism.lt.0) phism=phism+pi2
      if(phism.ge.pi2) phism=phism-pi2
      go to 2
    1 continue
      phid=phism+fi
      if(phid.ge.pi2) phid=phid-pi2
      if(phid.lt.0) phid=phid+pi2
    2 continue
c     print 900,ut,del,phid,phism,j
      return
      end
