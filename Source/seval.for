      real function seval(n,u,x,y,b,c,d)
      real x(n),y(n),b(n),c(n),d(n)
      data i/1/,k/1/
      if(i.ge.n)i=1
      if(u.lt.x(i)) go to 10
      if(u.le.x(i+1)) go to 30
   10 i=1
      j=n+1
   20 k=(i+j)/2
      if( u.lt.x(k) ) j=k
      if(u.ge.x(k)) i=k
      if(j.gt.(i+1)) go to 20
   30 dx=u-x(i)
      seval=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      k=k+1
      return
      end
