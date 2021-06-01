      subroutine ficon(an1,r,n,n1,n2,dt
     *                ,ddolg,vj,dtet,vi)

      dimension an1(n1,n2,n),vi(n1,n2,n),
     *          vj(n1,n2,n),r(n)
	allocatable rr(:)
      allocate (rr(n2))
      data pi/3.14159265359d0/,re/6.371e+8/
      rd=pi/180.
      it=n1-1
      do 2 i = 2 , it
        tet=(i-1)*dtet*rd
        ddr=ddolg*rd
        do 4 k=2,n
          do 3 j = 1 , n2
            j1=j-1
            j2=j+1
            if(j.eq.1) j1=n2
            if(j.eq.n2) j2=1
            s1=alog(an1(i,j,k))
            rk=r(k)+re
            am=(vj(i,j,k)-abs(vj(i,j,k)))/rk/sin(tet)/2./ddr*dt
            ap=(vj(i,j,k)+abs(vj(i,j,k)))/rk/sin(tet)/2./ddr*dt
            b2=(an1(i,j,k)-an1(i,j1,k))*ap/an1(i,j,k)
     *      +(an1(i,j2,k)-an1(i,j,k))*am/an1(i,j,k)
c           c1=(vj(i,j2,k)-vj(i,j1,k))/2./ddr*dt/rk/sin(tet)
c         aa=dt/(dtet*rd)/2./rk/sin(tet)
c     dd=aa*(vi(i+1,j,k)*sin(tet+dtet*rd)-
c    *   vi(i-1,j,k)*sin(tet-dtet*rd))
c         c1=c1*.5+dd*.5
           c1=0.
            rr(j)=s1-b2-c1
  3       continue
          do 5 j = 1 , n2
            an1(i,j,k)=exp(rr(j))
  5       continue
  4     continue
  2   continue
      deallocate (rr)
      return
      end

