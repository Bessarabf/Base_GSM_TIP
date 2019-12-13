      subroutine wwt (readfl,nfile1,kpar,dolm,ddolgt,nomsl,
     *       ntsl,nl,kdf,ldor,isp,md,par,pole,nr,mast)
      dimension ntsl(nl),kdf(20),par(nr),pole(ldor/4),mast(40)
      logical readfl
      nfile=nfile1
      if(nfile.eq.3)print 101,nfile,readfl
  101 format(' wwt :   nfile=',i3,'  readfl=',l8)
      if(nfile.ne.14) then 
        if(mast(13).eq.0) nfile=13  ! если дрейф не учитывать, f5
      end if 
      nob1=0
      nob2=0
!
      ns=0
      do i=1,nl
        ns=ns+ntsl(i)
      end do
! 
     if(nfile.ne.3) then 
        if(dolm.gt.0.) then 
          nn=dolm/ddolgt
          nob2=kpar*ns*nn
        end if 
     end if
      if(md.eq.0) then
        if(nomsl.ne.1) then 
          ns1=0
          k=nomsl-1
          do i=1,k
            ns1=ns1+ntsl(i)
          end do
          nob1=kpar*ns1
        end if 
        lpar=ntsl(nomsl)*kpar
      else  
        lpar=kpar*ns
      end if
      nob=nob1+nob2
      call inpout(readfl,nfile,kdf,ldor,
     *       isp,nob,lpar,par,pole,nr)
      return
      end

