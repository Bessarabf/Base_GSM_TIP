!         Считывание данных Miyahara (или других)
      subroutine botcalc_L(pgl,nh,its,ids,kpars,uts,dtets,ddolgs,
     *           rads,gkoor,pril,kpa,ntime)
     
     
      dimension pgl(kpars,nh,its,ids),gkoor(2,its,ids),rads(nh)
   
      dimension pril(kpa,its,ids,ntime),mkp(7)

      allocatable p(:,:,:),gir(:),tes(:)

      data mkp/1,2,3,7,10,11,12/
      allocate (p(kpars,its,ids),gir(ids),tes(its))

  900 format (a25)
  901 format (a1)
  902 format (F5.1,4F8.2)
  903 format (F5.1,3X,E9.2,8X,E9.2,5X,3F8.1)
  904 format (5E11.3)

!      расчет момента времени
!     данные записаны с интервалом ntime час
	nDT=24/ntime        ! time interval
!      l=(uts-dts)/3600./nDT+1   ! number point of massive PRIL
      l=(uts-0.1)/3600./nDT+1   ! number point of massive PRIL
      if(l.gt.ntime) then
		print*, 'massive pril exeeded in botcalc ',l,' uts=',uts
!	stop
		l=ntime
      end if

	 
!  longitude
      do i=1,ids
        gir(i)=ddolgs*(i-1)
      enddo

!  latitude 
      do i=1,its
        tes(i)=dtets*(i-1)
      enddo

C  для концентраций                
      CN=2.96e14+7.95e13+8.5e10 ! summury density
      A1=7.95E13/CN                
      A2=2.96E14/CN                
      A3=8.5E10/CN                 
      AMP=48.12E-24                
      
C**********************************************************
      do j = 1,ids
        do i = 1,its
	! 
	     
	     prr=pril(1,i,j,l) !
      
		 p(1,i,j)=prr*a1/amp ! O2
           p(2,i,j)=prr*a2/amp ! N2
           p(3,i,j)=prr*a3/amp ! O
           p(7,i,j)=pril(2,i,j,l)        ! T
	     p(10,i,j)=pril(3,i,j,l)       !
           p(11,i,j)=pril(4,i,j,l)       !
           p(12,i,j)=pril(5,i,j,l)       !
        enddo
      enddo
      do ikp=1,7 ! interpolation
       kpar=mkp(ikp) 
       call intpa(tes,its,dtets,gir,ids,ddolgs,gkoor,p,  ! interpolation to geom 
     *            kpar,pgl,kpars,nh)
      enddo
	
      call bonPGL1(pgl,kpars,nh,its,ids)
	
      
      call noznew( gir,tes,kpars,nh,its,ids,pgl) ! vector to geomag coor

	print*,'botcalc: end of interpolation'

      deallocate (p,gir,tes)
      return
      end
  

      subroutine noznew(gir,tes,kpars,nh,its,ids,pgl)
      dimension gir(ids),pgl(kpars,nh,its,ids),tes(its)
       i=1
        do j=1,its
          do l=1,ids
            t=tes(j)
            f=gir(l)
            r=g11t31(f,t)
            r=-r
            s=sin(r)
            c=cos(r)
            st=pgl(11,1,j,l)
            sf=pgl(12,1,j,l)
            p1=st*c-sf*s
            pgl(11,1,j,l)=p1
            p2=sf*c+st*s
            pgl(12,1,j,l)=p2
          enddo
        enddo
      return
      end

      subroutine intpa(tes,its,dtets,gir,ids,ddolgs,gkoor,p,
     *            kpar,pgl,kpars,nh)
      dimension tes(its),gir(ids),gkoor(2,its,ids),p(kpars,its,ids),
     *          pgl(kpars,nh,its,ids)
  900 format (2e12.3)
      do j=1,ids
        do i=1,its
          tet =gkoor(1,i,j)
          dolg=gkoor(2,i,j)
c         print *,' i,tet,dolg',i,tet,dolg
          call find(its,tet, tes,in)
          call find(ids,dolg,gir,jn)
          tetin=tes(in)
          dx=(tet-tetin)/dtets
          dolgjn=gir(jn)
c         print *,' tetin,dolgjn',tetin,dolgjn
          dy=(dolg-dolgjn)/ddolgs
          jn1=jn+1
          if(jn.eq.ids)jn1=1
c         print *,' in,in+1,jn,jn+1',in,in+1,jn,jn1
          p1=p(kpar,in,jn)
          p2=p(kpar,in+1,jn)
          p3=p(kpar,in,jn1)
          p4=p(kpar,in+1,jn1)
          f1=p1+(p2-p1)*dx
          f2=p3+(p4-p3)*dx
          ff=f1+(f2-f1)*dy
c         print *,' p1,p2,p3,p4,ff',p1,p2,p3,p4,ff
          pgl(kpar,1,i,j)=ff
        enddo
      enddo
      return
      end

      subroutine bonPGL1(pgl,kpars,nh,its,ids)
      dimension pgl(kpars,nh,its,ids)
     *         ,inp(4)        
	data inp/1,2,3,7/
      i2=its-1
      k=1

	do 3 i=1,3
          np=inp(i)
         
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
        do 4 j=1,ids
          pgl(np,k,1,j)=unp
          pgl(np,k,its,j)=usp
    4   continue
    
    3 continue
      return
      end


