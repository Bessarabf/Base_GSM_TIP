c  version GSM  15.04.14  
	subroutine globr(pgl,par,kpars,ddolgs,dtets,nh,
     *                 kdf,ldor,isp,pole,nr,its,ids,mass)
      dimension pgl(kpars,nh,its,ids),par(kpars,nh,its)
     *         ,kdf(20),pole(ldor/4)
     *         ,mass(30)
      logical readfl
c 772 format (' globr - whod ',i3,2e15.5,5i7)
c     print 772,kpars,ddolgs,dtets,nh,ldor,nr,its,ids
      readfl=.true.
      nfile=5
      npgl=kpars*nh*its*ids

      call globrw (nfile,readfl,pgl,pole,kpars,ldor,nh,
     *            kdf,isp,npgl,its,ids,mass)
              
c 777 format (' globr - wihod ',3e13.4)
c     print 777,pgl(1,1,1,1),pgl(2,1,1,1),pgl(3,1,1,1)
      return
      end
