      subroutine gsmtoham(qJGSM,qJHAM,nh,its,ids,nHAM,ddolgs,dtets)
      dimension qJHAM(ids,nHAM,its),qJGSM(its,ids,nh)
!
      allocatable tes(:),gir(:),pHAM(:,:),GSM(:,:)
     *           ,aMAGkoor(:,:,:)      
      allocate (tes(its),gir(ids),pHAM(ids,its),GSM(its,ids)
     *          ,aMAGkoor(2,its,ids))

!     GEOdetic longitudes
      do j=1,ids
        gir(j)=ddolgs*(j-1)
      enddo
!     GEOdetic latitudes  
      do i=1,its
        tes(i)=dtets*(i-1)
      enddo
	
      do i=1,its
	  tet=tes(i)
        do j=1,ids
	  dolg=gir(j)
          call ggmraw(0,dolg,tet,dolm,tetm)
          aMAGkoor(1,i,j)=tetm
          aMAGkoor(2,i,j)=dolm
        end do
      end do

      do k=1,nHAM
          
	    do i=1,its
	      do j=1,ids
		   GSM(i,j)=qJGSM(i,j,k)
	      end do
	    end do
	    call intpaGSM(tes,its,dtets,gir,ids,ddolgs,aMAGkoor,pHAM,GSM)
            do i=1,its
              do j=1,ids
	           qJHAM(j,k,i)=pHAM(j,i)
	      end do
	    end do
      end do
      deallocate (tes,gir,pHAM,GSM,aMAGkoor)
	
      return 
      end

      subroutine intpaGSM(tes,its,dtets,gir,ids,ddolgs,aMAGkoor,
     *                    pHAM,GSM)
      dimension tes(its),gir(ids),aMAGkoor(2,its,ids),pHAM(ids,its),
     *          GSM(its,ids)
  900 format (2e12.3)
      do j=1,ids
        do i=1,its
          tet =aMAGkoor(1,i,j)
          dolg=aMAGkoor(2,i,j)

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
          p1=GSM(in,jn)
          p2=GSM(in+1,jn)
          p3=GSM(in,jn1)
          p4=GSM(in+1,jn1)
          f1=p1+(p2-p1)*dx
          f2=p3+(p4-p3)*dx
          ff=f1+(f2-f1)*dy
c         print *,' p1,p2,p3,p4,ff',p1,p2,p3,p4,ff
          pHAM(j,i)=ff
        enddo
      enddo
      s np=0.
      s sp=0.
      i2=its-1
        do 1 j = 1 , ids
         snp=snp+pHAM(j,2)
         ssp=ssp+pHAM(j,i2)
   1    continue
c
        unp=snp/ids
        usp=ssp/ids
        do 4 j=1,ids
          pHAM(j,1)=unp
          pHAM(j,its)=usp
    4   continue
 
      return
      end
      
      subroutine vector_geo(nh,its,ids,vi,vj,ddolgs,dtets)
      dimension gir(ids),tes(its),vi(ids,nh,its),vj(ids,nh,its)
       do k=1,nh
        do i=2,its-1
          t=dtets*(i-1)
          do j=1,ids
            f=ddolgs*(j-1)
            r=g11t31(f,t)
            r=-r
            s=sin(r)
            c=cos(r)
            st=vi(j,k,i)
            sf=vj(j,k,i)
            p1=st*c-sf*s
            vi(j,k,i)=p1
            p2=sf*c+st*s
            vj(j,k,i)=p2
          end do
        end do
      end do
      return
      end
      
      subroutine bonvecham(vi,vj,n,n1,n2)
      dimension
     *         vi(n2,n,n1),vj(n2,n,n1)
      data pi/3.1415926/
      np=n1-1
      dfi=pi*2.0/n2
      dtet=pi/np
      n6=n2/2
      cosin=cos(dtet)
      do 1 k=1,n
       do 2 j=1,n6
        j1=j+n6
        fi=(j-1)*dfi
        fi180=fi+pi
c    . . . north pole
        as=sin(cosin*fi)
        ac=cos(cosin*fi)
        asin pi=sin(cosin*fi180)
        acos pi=cos(cosin*fi180)
        a1=vj(j,k,2)*ac+
     *     vi(j,k,2)*as
        a2=vj(j1,k,2)*acos pi+
     *     vi(j1,k,2)*asin pi
        a=0.5*(a1+a2)
        b1=vj(j,k,2)*as-vi(j,k,2)*ac
        b2=vj(j1,k,2)*asin pi-
     *     vi(j1,k,2)*acos pi
        b=(b1+b2)*0.5
        vi(j,k,1)=a*sin(fi)-b*cos(fi)
        vi(j1,k,1)=-vi(j,k,1)
        vj(j,k,1)=a*cos(fi)+b*sin(fi)
        vj(j1,k,1)=-vj(j,k,1)
c    . . . south pole
        a1=vj(j,k,np)*ac-
     *     vi(j,k,np)*as
        a2=vj(j1,k,np)*acos pi-
     *     vi(j1,k,np)*asin pi
        a=0.5*(a1+a2)
        b1=-vj(j,k,np)*as-vi(j,k,np)*ac
        b2=-vj(j1,k,np)*asin pi-
     *     vi(j1,k,np)*acos pi
        b=(b1+b2)*0.5
        vi(j,k,n1)=-a*sin(fi)-b*cos(fi)
        vi(j1,k,n1)=-vi(j,k,n1)
        vj(j,k,n1)=a*cos(fi)-b*sin(fi)
        vj(j1,k,n1)=-vj(j,k,n1)
   2   continue
   1  continue
      return
      end
