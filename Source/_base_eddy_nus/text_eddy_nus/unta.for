      function unta(l,n,x,y)
      dimension x(*),y(*)
      s=0.
      if(l.lt.n-1) then
        do i=l+1,n-1
          a=x(i+1)-x(i-1)
          b=y(i)
	    c=a*b
	    s=s+c
!	    write(8,'(i5,1p4g13.3)') i,a,b,c,s 
!          s=s+(x(i+1)-x(i-1))*y(i)
        end do
      end if
      s=s+(x(l+1)-x(l))*y(l)
      s=s+(x(n)-x(n-1))*y(n)
      unta=s*.5
      return
      end
