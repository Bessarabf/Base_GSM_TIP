      subroutine glqint(pgl1,par,pari,kpars,nh,its,ids,kdf,ldor,pole,
     *                  isp,ddolgs,dtets,sole,solen,rads,gkoor,
     *                  uts,delta,nse,gins,ins,mass,dts,ps,park,
     *                  nr,ni,kpart,ks,ddolgt,int,ntsl,nl,parj,parj1,
     *                  ntr,qom,qmax,iqo,mast,vir,E0,FAE)
      dimension pgl1(kpars,nh,its,ids),gins(ins,nh,its,ids)
     *         ,par(kpars,nh,its),pari(ins,nh,its)
     *         ,sole(nse),solen(nse),rads(nh),gkoor(2,its,ids)
     *         ,parj(nh,its,ids),parj1(nh,its),ps(10)
     *         ,mass(30),kdf(20),pole(ldor/4),ntsl(nl),park(ks)
     *         ,qom(nl),mast(40)
      dimension E0(its,ids),FAE(its,ids)
      dimension vir(nh,its,ids)
      logical readfl

      allocatable cO2plus(:,:,:),cNOplus(:,:,:)
      allocate (cO2plus(nh,its,ids),cNOplus(nh,its,ids))
      
      data key/0/

	readfl=.false.
      if(mass(20).eq.1) then
         call molio2(pgl1,gins,rads,vir,kpars,nh,its,dts,ids,
     *                ddolgs,ins,dtets,ntr)
      end if
  	print *,' glqint '
      if(mass(20).eq.2) then
      !
         open(77,file='molion.dat',form='unformatted')
         if(key.ne.0) then
           rewind(77)
	   read(77)cO2plus
           read(77)cNOplus  
	 end if
      end if    
      do 1 j = 1 , ids
       dolg=ddolgs*(j-1)
       
       DO i = 1 , its
        do k = 1 , nh
         do np = 1 , kpars
           par(np,k,i)=pgl1(np,k,i,j)
         end do
         do nin = 1 ,ins
          pari(nin,k,i)=gins(nin,k,i,j)
         end do
        end do
       END DO
       call ioniz_AE(sole,solen,rads,par,parj1,gkoor,uts,dtets,dolg,
     *               ddolgs,delta,nh,its,ids,nse,kpars,mass,ps,E0,FAE)
       if(mass(20).eq.0) then
         call molion (par,kpars,nh,its,pari,ins,mass,dts,ntr)
       else if(mass(20).eq.2) then
!	   call molio3S(cO2plus,cNOplus,par,ids,its,nh,kpars,pari,ins,
!     *                mass,dts,j,ntr,key)  
         call molio3nIT(cO2plus,cNOplus,par,ids,its,nh,kpars,pari,ins,
     *                  mass,dts,j,ntr,key)  
       end if
       call temol (par,pari,rads,mass,kpars,nh,its,dts,ntr,ins)
c . . . обход интерполяции шар-трубка при фиксировании ионосферы
      if(mass(13).ne.0) then
         call inst (dolg,ntsl,nl,ntr,ddolgt,kdf,ldor,isp,
     *              par,pari,pole,nr,ni,park,ks,int,rads,
     *              nh,its,dtets,kpars,qom,qmax,iqo,mast)
      end if
      nfile=9
      do  np = 1 , kpars
       do  k = 1 , nh
        do  i = 1 , its
         pgl1(np,k,i,j)=par(np,k,i)
        end do
       end do 
       end do

      do k = 1,nh
         do i = 1,its
           parj(k,i,j)=parj1(k,i)
         end do
	  end do
    1 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of longitude 

	
      key=1
      ! pole smoothing
      call botite(pgl1,kpars,nh,its,ids)
      if(mass(20).eq.2) then
!       call bonPGL_np(pgl1,kpars,nh,its,ids,6)
!		call bonPGL_np(pgl1,kpars,nh,its,ids,8)
                rewind(77)
		write(77) cO2plus
	        write(77) cNOplus
	        close(77)
        end if
      deallocate (cO2plus,cNOplus)
      return
      end

	subroutine bonPGL_np(pgl,kpars,nh,its,ids,np)
	! smoothing in polar for np- number parametr
      dimension pgl(kpars,nh,its,ids)
      i2=its-1
      do k=1,nh
         
c      ssp - sum s.pole
c      snp - sum n.pole
        s np=0.
        s sp=0.
        do 1 j = 1 , ids
         snp=snp+pgl(np,k,2,j)
         ssp=ssp+pgl(np,k,i2,j)
   1    continue
c
        unp=snp/ids
        usp=ssp/ids
        do  j=1,ids
          pgl(np,k,1,j)=unp
          pgl(np,k,its,j)=usp
        end do 
    
      end do
      return
      end
