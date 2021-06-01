c . . . ver. 16.04.14 - add idt to interface parameters
      subroutine readinf(pole,ldor,god,day,ut0,ut1,utk,dtt,dts,tau,kdf,
     *                   kdu,kpart,kpars,int,ins,its,ids,idt,
     *                   ddolgt,dtett,ddolgs,dtets,ntr,nv,nh,nl,
     *                   ntsl,ntpl,ns,fa0,fs,nsu,nse,
     *                   bmpy,bmpz,vsol,csol,ap0,pkp0,dst0,
     *                   ae0,al0,au0,gamma,rmaxt,b,c,dh,rmin,
     *                   tet1,tet2,tet3,eps0,om,
     *                   fac1,fac2,fac3,pdpc,
     *                   nzapt,nzaps,nadrt,nadrs)
      dimension kdf(20),kdu(20),ntsl(*),pole(1024)
      integer god,day
c        reading inform area
c
      ldor=  pole(1) 
      god =  pole(2) 
      day =  pole(3) 
      ut0 =  pole(4) 
      ut1 =  pole(5) 
      utk =  pole(6) 
      dtt =  pole(7) 
      dts =  pole(8) 
      tau =  pole(9) 
      do i=1,20
        kdf(i)=pole(i+9) 
        kdu(i)= pole(i+29)
      end do

!     число записей итерационных области шара и трубки 
!     должны совпадать с соответстующими областями шара и трубки
	if (kdu(11).ne.kdu(5)) then
	    kdu(11)=kdu(5)
	    kdu(12)=kdu(6)
	    kdu(13)=kdu(6)
	    kdu(14)=kdu(6)
	    do kk=12,15
	      kdf(kk)=kdu(kk-1)+kdf(kk-1)
	    end do
	end if
      kpart =pole(50) 
      kpars =pole(51) 
      int   =pole(52) 
      ins   =pole(53) 
      ddolgt=pole(54) 
      dtett =pole(55) 
      ddolgs=pole(56)
       
      dtets =pole(57)
      its=180.1/dtets+1
      ids=360.1/ddolgs 
	idt=360.1/ddolgt 
      ntr   =pole(58) 
      nv    =pole(59) 
      nh    =pole(60) 
      nl    =pole(61) 
      ntpl=0
      do i=1,nl
        ntsl(i)=pole(i+61)
        ntpl=ntpl+ntsl(i)
      end do
      
      ns = pole(226)
      fa0= pole(227)
      fs = pole(228)
      nsu= pole(229)
      nse= pole(230)
      
      gamma= pole(385)
      rmaxt= pole(386)
      b    = pole(387)
      c    = pole(388)
      dh   = pole(389)*1.e5
      rmin = pole(390)
      tet1 = pole(391)
      tet2 = pole(392)
      tet3 = pole(393)
      eps0 = pole(394)
      om   = pole(395)
      fac1 = pole(396)
      fac2 = pole(397)
      fac3 = pole(398)
      pdpc = pole(399)

      nzaps=  pole(600)
      nadrs=  pole(601)
      nzapt=  pole(602)
      nadrt=  pole(603)
      return
      end
