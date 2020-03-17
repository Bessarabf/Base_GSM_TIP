      subroutine eiscat(etet,edol,eheig,namefl,par1,par2,potef,pole,
     *                  kpars,kpart,nh,its,ntr,
     *                  idt,nl2,mdor,ntsl,nl,kdf,kdu,isp,tau,ks,
     *                  ddolg,dtets,dtet,nzap,adres,nr,ut,park,dts,
     *                  nrec,ip,nrecl,mast)

      parameter(kpars16=16)
! kpars16 - only 16 parameters
      character*25 namefl
      logical readfl/.true./
      integer ntsl(nl),kdf(20),adres,kdu(20)
      
      real hm(70),cm1(70,2,8)
     *    ,fil(8,70)
      
      integer npar(8)/1,2,3,4,5,6,7,8/
      dimension par1(kpars,nh,its), par2(kpars,nh,its)
     *         ,pole(mdor/4), potef(ntr,idt,nl2),park(ks)
     *         ,nf(17),mast(40)
     
      allocatable a(:),cmi(:,:),ssh(:), out(:,:)

	
      common /flstat/ nf
      data  imin/1/, ns/1/
      data md/1/
	
      allocate (a(its),cmi(its,nh),ssh(its),out(its,nh))

      ntpl=ks/2
!
! ifix - real into intefer
!
      dolg1=ifix(edol/ddolg)*ddolg
      if(dolg1.eq.360.) dolg1=0.
      dolg2=dolg1+ddolg
      if(dolg2.eq.360.) dolg2=0.
      open (12,file=namefl,access='direct',recl=nrecl)
      nfile = 5 ! record numder in kdf() 
! first plosk
      call wws (readfl,nfile,kpars,dolg1,ddolg,
     *          tet,dtets,nh,kdf,mdor,isp,md,par1,pole,nr)
! second plosk
      call wws (readfl,nfile,kpars,dolg2,ddolg,
     *          tet,dtets,nh,kdf,mdor,isp,md,par2,pole,nr)
      coef2=(edol-dolg1)/ddolg
      coef1=1.-coef2
      do  k=1,its
        do j=1,nh
         do i=1,6
           if(par1(i,j,k).lt.1.e-3) par1(i,j,k)=1.e-3
           if(par2(i,j,k).lt.1.e-3) par2(i,j,k)=1.e-3
           par1(i,j,k)=alog10(par1(i,j,k))
           par2(i,j,k)=alog10(par2(i,j,k))
         end do 
         do i=13,16
            if(par1(i,j,k).lt.1.e-3) par1(i,j,k)=1.e-3
            if(par2(i,j,k).lt.1.e-3) par2(i,j,k)=1.e-3
            par1(i,j,k)=alog10(par1(i,j,k))
            par2(i,j,k)=alog10(par2(i,j,k))
         end do 
! linear interpolation
         do i=1,kpars16
            par1(i,j,k)=par1(i,j,k)*coef1+par2(i,j,k)*coef2
         end do
        end do 
      end do

     
      do 2 i=1,its
            a(i)=90.-(i-1)*dtets
        if ( etet.ge.a(i) ) then
          i1=i
          i2=i-1
          coef2=( etet-a(i1) )/dtets
          coef1=1.-coef2
          goto 3
        endif
 2    continue
 3    continue
      do 4 n=1,nh
         do 5 k=1,kpars16
            out(k,n)=par1(k,n,i1)*coef1+par1(k,n,i2)*coef2
            if(k.le.6.or.k.ge.13) out(k,n)=10.**(out(k,n))
 5       continue
 4    continue
      d0=edol
      tet=90.-etet
      idm=d0/ddolg
      coef1=(d0-idm*ddolg)/ddolg
      coef2=1.-coef1
      do 7 i=1,2
         id=idm+i-1
         if(id.ge.idt)id=id-idt
         dolg=ddolg*id
         nfile = 6
         call wwt(.true.,nfile,kpart,dolg,ddolg,i,ntsl,nl,kdf,mdor,
     *            isp,1,par1,pole,nr,mast)
        do i1=1,kpart
          call plosk(par1,park ,hm,cmi,ssh,tet,cm1(1,i,i1),ns,
     *               ip,imin,ntsl,ntpl,nl,npar(i1),its,nh,nl2)
 7      end do
      end do 
      do  i1=1,kpart
         do  id=1,ip
           fil(i1,id)=cm1(id,1,i1)*coef2+cm1(id,2,i1)*coef1
         end do 
      end do 
      do j=1,ip
        do i=1,3
           fil(i,j)=10**fil(i,j)
        end do
      end do
! potential
        nfile = nf(15)
        n3=ntr*idt*nl2
        call wpotef(.true.,potef,n3,kdf,kdu,mdor,isp)
        id1=idm+1
        if(id1.eq.idt+1)id1=1
        idd1=id1+1
        if(idd1.gt.idt)idd1=idd1-idt
        it1=tet/dtet+1
          rip = ip
        write(12,rec=nrec) ut,etet,edol,eheig,rip,
     *  ((OUT(I,J),I=1,KPARS16),J=1,NH),((FIL(I,J),I=1,KPART),J=1,Ip),
     *  potef(ntr,id1,it1),potef(ntr,id1,it1+1),potef(ntr,idd1,it1)
     *  ,potef(ntr,idd1,it1+1)
        close (12)
!        write(*,'(a8,a21,a8)')' end of ',namefl,' filling'
	deallocate (a,cmi,ssh,out)
      return
      end

! subrout for eiscat
!
      subroutine plosk(pk,pk3,hm,cmi,ssh,se,cm1,ns,ip,imin,ntsl
     *                ,ntpl,nl,npar,its,nh,nl2)

      integer ntsl(nl)
      real hm(ip),cm1(ip,ns),pk(8,ntpl),
     *     pk3(2,ntpl),cmi(its,nh),ssh(its),se(ns)
      allocatable sm(:),cm(:)
      allocate (sm(nl2-2),cm(nl2-2))      
      ih=imin
      do 29 ihp=1,ip
        hm(ihp)=pk3(1,ih)*1.e-5
        i0=0
        i=0
        do jt=1,nl
           n2=ntsl(jt)/2
	   if(ih.gt.n2)goto 19
	   i=i+2
	   ih1=i0+ih
	   sm(i)=pk3(2,ih1)
	   cm(i)=cf(pk,ntpl,8,npar,ih1)
	   i0=i0+ntsl(jt)
	   ih1=i0+1-ih
	   sm(i-1)=pk3(2,ih1)
	   cm(i-1)=cf(pk,ntpl,8,npar,ih1)
        end do
  19    call turn1(sm,cm,i)
c               call spline(i,sm,cm,b,c,d)
        do j=1,ns
          cm1(ihp,j)=val(i,se(j),sm,cm)
          if(ih.le.15.and.npar.eq.9) cm1(ihp,j)=
     *                               alog10(10.**(cm1(ihp,j))+
     *                    10.**val(19,se(j),ssh,cmi(1,ih+15)))
        end do 
  29  ih=ih+1
      deallocate (sm,cm)
      return
      end
c
      function cf(pk,ntpl,npk,npar,j1)
      real pk(npk,ntpl)
      if(npar.ne.9)then
        cf=pk(npar,j1)
        if(npar.ge.4.and.npar.le.6) cf=cf*pk(npar-3,j1)
        if(npar.le.3)then
c       if(npar.le.3.or.npar.ge.7)then
           cf=alog10(cf)
cc         if(cf.lt.0.) cf=0.
           if(cf.lt.-3.) cf=-3.
        endif
      else
        cf=alog10(pk(1,j1)+pk(2,j1))
      endif
      return
      end
c
      subroutine turn1(h,c,n)
      real h(n),c(n)
      j=1
  3   if(h(j).le.h(j+1))goto 2
      call turn(h(j))
      call turn(c(j))
      j=j-2
      if(j.eq.-1)j=0
  2   j=j+1
      if(j.lt.n)goto 3
      return
      end
c
      subroutine turn(sm)
      real sm(2)
      a=sm(1)
      sm(1)=sm(2)
      sm(2)=a
      return
      end
cc
      real function val(n,u,x,y)
      real x(n),y(n)
      data i/1/,k/1/
      if(i.ge.n)i=1
      if(u.lt.x(i)) go to 10
      if(u.le.x(i+1)) go to 30
   10 i=1
      j=n+1
   20 k=(i+j)/2
      if( u.lt.x(k) ) j=k
      if(u.ge.x(k)) i=k
      if(j.gt.(i+1)) go to 20
   30 dx=(u-x(i))/(x(i+1)-x(i))
      val=y(i)+dx*(y(i+1)-y(i))
      return
      end
