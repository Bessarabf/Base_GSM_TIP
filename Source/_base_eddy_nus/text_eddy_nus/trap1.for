      subroutine trap1(f,x,n,q,ki,ke)
c     . . . quadrating programm from x(ki) to x(ke)
      dimension f(n),x(n)
      q=0.
      do 1 i=ki,ke-1
       q=q+(f(i)+f(i+1))*(x(i+1)-x(i ))*.5
   1  continue
      return
      end


