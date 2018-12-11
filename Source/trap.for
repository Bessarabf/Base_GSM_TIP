      subroutine trap(f,x,n,stra)
      dimension f(n),x(n)
      n 1=n-1
      stra=0.
      do 1 k=1,n 1
       chisl=f(k)+f(k+1)
       znam=x(k+1)-x(k)
       stra=stra+0.5*chisl*znam
   1  continue
      return
      end
