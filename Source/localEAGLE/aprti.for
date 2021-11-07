      subroutine aprti(nx,tn,ti)
      
      dimension tn(*),ti(*)
      do 1 i=1,nx
c       ti(i)=tn(i)*5.
        ti(i)=tn(i)
    1 continue
      return
      end
