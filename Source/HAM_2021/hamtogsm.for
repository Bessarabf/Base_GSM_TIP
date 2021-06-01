      subroutine hamtogsm(tHAM,gsmHAM,gkoor,nh,its,ids,nHAM,	
     *                    ddolgs,dtets)
      allocatable tes(:),gir(:),pHAM(:,:),GSM(:,:)
      dimension tHAM(ids,nHAM,its),gsmHAM(nh,its,ids),gkoor(2,its,ids)
      allocate (tes(its),gir(ids),pHAM(ids,its),GSM(its,ids))
	do i=1,ids
        gir(i)=ddolgs*(i-1)
      enddo
!  latitudes massive 
      do i=1,its
        tes(i)=dtets*(i-1)
      enddo
      do k=1,nHAM
          kh=nHAM-k+1 ! reverse order
	    do i=1,its
	      do j=1,ids
		   pHAM(j,i)=tHAM(j,kh,i)
	      end do
	     end do
	     call intpaHAM(tes,its,dtets,gir,ids,ddolgs,gkoor,pHAM,GSM)
          do i=1,its
            do j=1,ids
		     gsmHAM(k,i,j)=GSM(i,j)
	      end do
	    end do
      end do
      deallocate (tes,gir,pHAM,GSM)
      return 
      end

      subroutine intpaHAM(tes,its,dtets,gir,ids,ddolgs,gkoor,pHAM,GSM)
      dimension tes(its),gir(ids),gkoor(2,its,ids),pHAM(ids,its),
     *         GSM(its,ids)
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
c         print *,' tetin,dolgjn',tetin, 
          dy=(dolg-dolgjn)/ddolgs
          jn1=jn+1
          if(jn.eq.ids)jn1=1
c         print *,' in,in+1,jn,jn+1',in,in+1,jn,jn1
          p1=pHAM(jn,in)
          p2=pHAM(jn,in+1)
          p3=pHAM(jn1,in)
          p4=pHAM(jn1,in+1)
          f1=p1+(p2-p1)*dx
          f2=p3+(p4-p3)*dx
          ff=f1+(f2-f1)*dy
c         print *,' p1,p2,p3,p4,ff',p1,p2,p3,p4,ff
          GSM(i,j)=ff
        enddo
      enddo
      s np=0.
      s sp=0.
      i2=its-1
        do 1 j = 1 , ids
         snp=snp+GSM(2,j)
         ssp=ssp+GSM(i2,j)
   1    continue
c
        unp=snp/ids
        usp=ssp/ids
        do 4 j=1,ids
          GSM(1,j)=unp
          GSM(its,j)=usp
    4   continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine vecGSM(UgsmHAM,VgsmHAM,nh,its,ids,nHAM,ddolgs,dtets) 
      
      dimension UgsmHAM(nh,its,ids),VgsmHAM(nh,its,ids)
       
       do k=1,nham
        t=0.    ! colatitude
        do j=1,its
          f=0.
          do l=1,ids
            !t=tes(j)
            !f=gir(l)
            r=g11t31(f,t)
            r=-r
            s=sin(r)
            c=cos(r)
            st=-UgsmHAM(k,j,l) ! revert meridional wind
            sf=VgsmHam(k,j,l)  ! zonal wind 
            p1=st*c-sf*s
            UgsmHAM(k,j,l)=p1*100. ! cm/s
            p2=sf*c+st*s
            VgsmHam(k,j,l)=p2*100.  ! cm/s
            f=f+ddolgs

          enddo
          t=t+dtets
        enddo
      end do
      return
        end
        