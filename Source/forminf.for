      subroutine forminf(ldor,god,day,ut0,ut1,utk,dtt,dts,tau,kdf,
     *                   kdu,kpart,kpars,int,ins,
     *                   ddolgt,dtett,ddolgs,dtets,ntr,nv,nh,nl,
     *                   ntsl,rads,mass,mast,mas,ns,fa0,fs,nsu,nse,
     *                   sole,solet,solen,solu,
     *                   bmpy,bmpz,vsol,csol,ap0,pkp0,dst0,
     *                   ae0,al0,au0,gamma,rmaxt,b,c,dh,rmin,
     *                   tet1,tet2,tet3,eps0,om,
     *                   fac1,fac2,fac3,pdpc,pril,
     *                   nzapt,nzaps,nadrt,nadrs)
      dimension kdf(20),kdu(20),ntsl(44),pole(1024),pril(*)
      dimension sole(40),solet(40),solen(40),solu(24)
      dimension mas(10),mass(30),mast(40),rads(30)
      integer god,day
      character*8 namemes,im*80
  900 format(a25)
  901 format(3(5i6/),4(10i6/),2i6,2e10.3)
c        Формирование информационной области
c
      pole(1)  = ldor
      pole(2)  = god
      pole(3)  = day
      pole(4)  = ut0
      ut1=ut0
      pole(5)  = ut1
      pole(6)  = utk
      pole(7)  = dtt
      pole(8)  = dts
      pole(9)  = tau
      do i=1,20
        pole(i+9) =kdf(i)
        pole(i+29)=kdu(i)
      end do
      pole(50) = kpart
      pole(51) = kpars
      pole(52) = int
      pole(53) = ins
      pole(54) = ddolgt
      pole(55) = dtett
      pole(56) = ddolgs
      pole(57) = dtets
      pole(58) = ntr
      pole(59) = nv
      pole(60) = nh
      pole(61) = nl
      do i=1,nl
        pole(i+61)=ntsl(i)
      end do
      do i=1,nh
        pole(i+105)=rads(i)
      end do
      do i=1,30
        pole(i+135)=mass(i)
      end do
      do i=1,40
        pole(i+165)=mast(i)
      end do
      do i=1,10
        pole(i+205)=mas(i)
      end do
      pole(226)= ns
      pole(227)= fa0
      pole(228)= fs
      pole(229)= nsu
      pole(230)= nse
      do i=1,nse
        pole(i+230) = sole(i)
        pole(i+270) = solet(i)
        pole(i+311) = solen(i)
      end do
      do i=1,nsu
        pole(i+350) = solu(i)
      end do
      pole(375)= bmpy
      pole(376)= bmpz
      pole(377)= vsol
      pole(378)= csol
      pole(379)= ap0
      pole(380)= pkp0
      pole(381)= dst0
      pole(382)= ae0
      pole(383)= al0
      pole(384)= au0
      pole(385)= gamma
      pole(386)= rmaxt
      pole(387)= b
      pole(388)= c
      pole(389)= dh*1.e-5
      pole(390)= rmin
      pole(391)= tet1
      pole(392)= tet2
      pole(393)= tet3
      pole(394)= eps0
      pole(395)= om
      pole(396)= fac1
      pole(397)= fac2
      pole(398)= fac3
      pole(399)= pdpc
      do i=1,143
        pole(i+399) = pril(i)
      end do
      sum=0
      pole(600)= nzaps
      pole(601)= nadrs
      pole(602)= nzapt
      pole(603)= nadrt
      do i=1,20
      sum = sum+kdu(i)
      end do
      pole(1000) =sum
      pole(1004) = 4
      pole(1005) = 5
      pole(1006) = 6
      pole(1007) = 7
      pole(1008) = 8
      pole(1009) = 9
      pole(1010) = 10
      pole(1011) = 11
      pole(1012) = 12
      pole(1013) = 13
      pole(1014) = 14
      write(4,rec=1)pole
      call calcdat(ut0,ut1,day,god,ng,nden,namemes,np,min)
  920 format(i19,' year',i3,1x,a8,i3,' hour.',i2,' min.'/)
      print   920,ng,nden,namemes,np,min
	print*,' god,day ', god,day

      return
      end
