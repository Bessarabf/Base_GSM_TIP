      subroutine spline(n,x,y,b,c,d)
      real x(n),y(n),b(n),c(n),d(n)
      data k/1/
      nm1=n-1
      if(n.lt.2) go to 999
      if(n.lt.3) go to 50
      d(1)=x(2)-x(1)
      c(2)=(y(2)-y(1))/d(1)
      do 10 i=2,nm1
        d(i)=x(i+1)-x(i)
        b(i)=2.*(d(i-1)+d(i))
        c(i+1)=(y(i+1)-y(i))/d(i)
        c(i)=c(i+1)-c(i)
   10 continue
      b(1)=-d(1)
      b(n)=-d(n-1)
      c(1)=0.
      c(n)=0.
      if(n.eq.3) go to 15
      c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1)=c(1)*d(1)**2/(x(4)-x(1))
      c(n)=-c(n)*d(n-1)**2/(x(n)-x(n-3))
   15 continue
      do 20 i=2,n
        t=d(i-1)/b(i-1)
        b(i)=b(i)-t*d(i-1)
        c(i)=c(i)-t*c(i-1)
   20 continue
      c(n)=c(n)/b(n)
      do 30 ib=1,nm1
        i=n-ib
        c(i)=(c(i)-d(i)*c(i+1))/b(i)
   30 continue
      b(n)=(y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.*c(n))
      do 40 i=1,nm1
        b(i)=(y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.*c(i))
        d(i)=(c(i+1)-c(i))/d(i)
        c(i)=3.*c(i)
   40 continue
      c(n)=3.*c(n)
      d(n)=d(n-1)
      go to 999
   50 continue
      b(1)=(y(2)-y(1))/(x(2)-x(1))
      c(1)=0.
      d(1)=0.
      b(2)=b(1)
      c(2)=0.
      d(2)=0.
  999 continue
      return
      end

