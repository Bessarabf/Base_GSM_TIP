      subroutine wind0(pgl,kpars,nh,its,ids)
c     . . . нулевые ветра на н.границе
      dimension pgl(kpars,nh,its,ids)
      do i=1,its
        do j=1,ids
         pgl(11,1,i,j)=0.
         pgl(12,1,i,j)=0.
        end do
       end do
       return
      end
                                               
