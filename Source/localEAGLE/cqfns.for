      subroutine cqfns(cqm,fa,hsl,jz,nol,nu,nvn,qs3,tsl,qfn,qfs)
      integer jz,nol,nu(*)
      real hsl(nvn,*),qs3(*),tsl(nvn,*)
      integer i,j,nj,ns1,nsn
      real a,c,ca,pd,r,re,sq3,tet,x(19),y(19)
      data re/6371.2/,pd/1.74532925e-2/,sq3/1.7320508/
      data c0/-201.69156/,c1/12.4424/,c2/-0.10053425/,c3/0.000352002/
      data c4/-4.32964e-7/
	ns1=1
      r=re+hsl(ns1,nol)
      a=tsl(ns1,nol)*pd
      s=sin(a)
      scx=r/(s*s*2.)*1e5
      h=(((c4*fa+c3)*fa+c2)*fa+c1)*fa+c0
      k=0
      do i=3,jz
        if(h.lt.hsl(i,nol)) then
          k=i
          exit
        end if
      end do
      a=hsl(k-1,nol)
      s=(h-a)/(hsl(k,nol)-a)
      t=s*(tsl(k,nol)-tsl(k-1,nol))+tsl(k-1,nol)
      i=1
      tet=t*pd
      c=cos(tet)
      ca=sqrt(c*c*3.+1.)
      a=c*sq3+ca
      x(i)=(alog(a)/sq3+ca*c)*scx
      y(i)=s*(qs3(k+1)-qs3(k))+qs3(k)
      do j=k,jz
        i=i+1
        tet=tsl(j,nol)*pd
        c=cos(tet)
        ca=sqrt(c*c*3.+1.)
        a=c*sq3+ca
        x(i)=(alog(a)/sq3+ca*c)*scx
        y(i)=qs3(j)
      end do
      qfs=unta(1,i,x,y)*cqm+1.e7
      j=1
      nj=nu(nol)
      y(j)=s*(qs3(nj-k+1)-qs3(nj-k+2))+qs3(nj-k+2)
      nsn=nj-jz+1
      do i=nj-k+1,nsn,-1
        j=j+1
        y(j)=qs3(i)
      end do
      qfn=unta(1,j,x,y)*cqm+1.e7
      return
      end