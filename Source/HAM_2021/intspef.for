      subroutine intspef(idol,par,nr,rads,nh,park,
     *        ks,ntsl,nl,
     *        dtett,u,nui,ddolgs,dtets,pot,na,ids,ncn,nc)
      dimension ntsl(nl),par(nr),rads(nh),park(ks),pot(na,ids,ncn),
     *        msum(45),grz(16),u(nui)
      data pi/3.1415926/,re/6371.02e5/
      cr=180./pi
      msum(1)=0
      do nomsl=2,nl
        nm=nomsl-1
        msum(nomsl)=msum(nm)+ntsl(nm)
      end do
      j=nl
      i1=msum(j)+ntsl(j)/2
      x1=rads(1)
      x2=park(i1*2+1)
      i2=msum(j)+ntsl(j)
      i3=(ncn+1)/2
      vn1=par(i2+1)
      vn2=par(i1+1)
      nm=na-1
      do i=2,nm
        xi=rads(i)
c        if(xi.gt.x2)then
	  do while (xi.gt.x2)
          j=j-1
          x1=x2
          i1=msum(j)+ntsl(j)/2
          x2=park(i1*2+1)
          vn1=vn2
          vn2=par(i1+1)
c        end if
        end do
        call int2pef(x1,x2,xi,vn1,vn2,plm)
        pot(i,idol,i3)=plm
      end do
      k1=ncn/2
      ncn1=(ncn+1)/2
      nui1=nui-1
      do i=1,nm
        j=2
        i1=msum(1)+i
        i2=msum(j)+i
        x2=park(i1*2)
        x1=park(i2*2)
        vn2=par(i1)
        vn1=par(i2)
        do kk=2,k1
          k=ncn-kk
          xi=dtett*k
c          if(xi.lt.x1)then
          do while (xi.le.x1)	
            j=j+1
            i1=i2
            x2=x1
            vn2=vn1
            if(j.gt.nui1)then
              x1=90.
              vn1=pot(i,idol,ncn1)
            else
              i3=msum(j)+(ntsl(j)+1)/2
              if(park(i1*2-1).gt.park(i3*2-1))then
                x1=90.
                vn1=pot(i,idol,ncn1)
              else
                i2=msum(j)+i
                x1=park(i2*2)
                vn1=par(i2)
              end if
            end if
c          end if
	    end do
c          if(xi.gt.x2)then
	    do while (xi.gt.x2) 
            i2=i1
            j=j-1
            i1=msum(j)+i
            x1=x2
            x2=park(i1*2)
            vn1=vn2
            vn2=par(i1)
c          end if
	    end do	 
          call int2pef(x1,x2,xi,vn1,vn2,plm)
          pot(i,idol,k+1)=plm
        end do
      end do
      do i=1,nm
        do k=5,18
          kk=ncn-k+1
          pot(i,idol,k)=pot(i,idol,kk)
        end do
      end do
      do i=1,nm
        j=2
        i1=msum(1)+ntsl(1)-i+1
        i2=msum(2)+ntsl(2)-i+1
        x2=park(i2*2)
        x1=park(i1*2)
        vn2=par(i2)
        vn1=par(i1)
        do kk=2,4
          xi=dtett*(kk-1)
c          if(xi.lt.x1)then
          do while(xi.le.x1)	
            j=j-1
            i2=i1
            x2=x1
            vn2=vn1
            i1=msum(j)+ntsl(j)-i+1
            x1=park(i1*2)
            vn1=par(i1)
c          end if
c          if(xi.gt.x2)then
	    end do
          do while(xi.gt.x2)
            i1=i2
            j=j+1
            i2=msum(j)+ntsl(j)-i+1
            x1=x2
            x2=park(i2*2)
            vn1=vn2
            vn2=par(i2)
c          end if
	    end do
          call int2pef(x1,x2,xi,vn1,vn2,plm)
          pot(i,idol,kk)=plm
        end do
      end do
      return
      end