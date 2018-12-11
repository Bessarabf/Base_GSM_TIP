      subroutine aprhe(nx,cihe,vihe)
      dimension cihe(*),vihe(*)
      do 1 i=1,nx
        cihe(i)=1.e-3
        vihe(i)=0.
    1 continue
      return
      end
