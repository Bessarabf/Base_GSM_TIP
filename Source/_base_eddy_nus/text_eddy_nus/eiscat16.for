      subroutine eiscat(etet,edol,eheig,namefl,par1,par2,potef,pole,
     *                  kpars,kpart,nh,its,ntr,
     *                  idt,nl2,mdor,ntsl,nl,kdf,kdu,isp,tau,ks,
     *                  ddolg,dtets,dtet,nzap,adres,nr,ut,park,dts,
     *                  nrec,ip,nrecl,mast)

	parameter(kpars16=16)
	!kpars16 - only 16 parameters
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
      dolg1=ifix(edol/ddolg)*ddolg
      if(dolg1.eq.360.)dolg1=0.
      dolg2=dolg1+ddolg
      if(dolg2.eq.360.)dolg2=0.
      open (12,file=namefl,access='direct',recl=nrecl)
          nfile = 5
      call wws (readfl,nfile,kpars,dolg1,ddolg,
     *          tet,dtets,nh,kdf,mdor,isp,md,par1,pole,nr)
      call wws (readfl,nfile,kpars,dolg2,ddolg,
     *          tet,dtets,nh,kdf,mdor,isp,md,par2,pole,nr)
          coef2=(edol-dolg1)/ddolg
          coef1=1.-coef2
      do 1 k=1,its
        do 1 j=1,nh
         do 10 i=1,6
          if(par1(i,j,k).lt.1.e-3) par1(i,j,k)=1.e-3
          if(par2(i,j,k).lt.1.e-3) par2(i,j,k)=1.e-3
          par1(i,j,k)=alog10(par1(i,j,k))
          par2(i,j,k)=alog10(par2(i,j,k))
 10      continue
         do 11 i=13,16
          if(par1(i,j,k).lt.1.e-3) par1(i,j,k)=1.e-3
          if(par2(i,j,k).lt.1.e-3) par2(i,j,k)=1.e-3
          par1(i,j,k)=alog10(par1(i,j,k))
          par2(i,j,k)=alog10(par2(i,j,k))
 11      continue
         do 1 i=1,kpars16
            par1(i,j,k)=par1(i,j,k)*coef1+par2(i,j,k)*coef2
 1      continue
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
 5        continue
 4      continue
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
     *           isp,1,
     *           par1,pole,nr,mast)
        do 7 i1=1,kpart
!          call plosk(par1,park ,hm,cmi,ssh,tet,cm1(1,i,i1),ns,
!     *               ip,imin,ntsl,ntpl,nl,npar(i1),its,nh,nl2)
 7      continue
      do 8 i1=1,kpart
        do 8 id=1,ip
          fil(i1,id)=cm1(id,1,i1)*coef2+cm1(id,2,i1)*coef1
 8    continue
          do 9 j=1,ip
          do 9 i=1,3
            fil(i,j)=10**fil(i,j)
 9        continue
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
