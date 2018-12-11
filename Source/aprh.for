      subroutine aprh(nx,cih,vih)
      dimension cih(*),vih(*)
      do 1 i=1,nx
        cih(i)=1.e-3
        vih(i)=0.
    1 continue
      return
      end
