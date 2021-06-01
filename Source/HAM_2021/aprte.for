      subroutine aprte(nx,tn,te)
      dimension tn(*),te(*)
      do 1 i=1,nx
c       te(i)=tn(i)*5.
        te(i)=tn(i)
    1 continue
      return
      end
