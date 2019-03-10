      subroutine find(n,u,x,m)
      dimension x(n)
      data i/1/
c     if(u.gt.x(n).or.u.le. x(1)) print 900,u,x(1),x(n)
c 900 format(' find'/
c    *' значение u выходит за пределы массива x в ',a4/
c    *'   u=',e10.3,'  x(1)=',e10.3,'  x(n)=',e10.3)
      if(i.ge.n) i=1
      if(u.lt.x(i)) go to 10
      if(u.le.x(i+1)) go to 30
c
   10 i=1
      j=n+1
   20 k=(i+j)/2
      if(u.lt.x(k)) j=k
      if(u.ge.x(k)) i=k
      if(j.gt.i+1) go to 20
   30 m=i
      return
      end
