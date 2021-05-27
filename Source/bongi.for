      subroutine bongi(gins,ins,nh,its,ids)
      dimension gins(ins,nh,its,ids)
      i2 = its - 1
      do i = 1, 3
        np = i
        if(i == 2) np = 5
        if(i == 3) np = 6
        do  k = 1 , nh
c      ssp - sum s.pole
c      snp - sum n.pole
          snp = 0.
          ssp = 0.
          do j = 1, ids
            snp = snp + gins(np,k,2,j)
            ssp = ssp + gins(np,k,i2,j)
          enddo ! j
          unp = snp / ids
          usp = ssp / ids
          do j = 1, ids
            gins(np,k,1,j) = unp
            gins(np,k,its,j) = usp
          enddo ! j
        enddo ! k
      enddo ! i
      end
