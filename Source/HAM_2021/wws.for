      subroutine wws (readfl,nfile,kpar,dolg,ddolgs,
     *       tet,dtets,nh,kdf,ldor,isp,md,par,pole,nr)
      dimension kdf(20),par(nr), pole(ldor/4)
      logical readfl
c     print 100
  100 format(' wws')
c     call timen
      nn=180./dtets+1
      lpar=nh*kpar
      nob1=0
      nob2=0
      if(dolg.ne.0.)nob2=lpar*nn*(dolg/ddolgs)
!
      if(md.eq.0) then 
        if(tet .ne.0.)nob1=(tet/dtets)*lpar
      else 
        lpar=lpar*nn
      end if
      nob=nob1+nob2
      call inpout(readfl,nfile,kdf,ldor,
     *       isp,nob,lpar,par,pole,nr)
c     call timen
      return
      end
