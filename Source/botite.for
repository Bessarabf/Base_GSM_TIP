      subroutine botite(pgl,kpars,nh,its,ids)
      dimension pgl(kpars,nh,its,ids)
      i2=its-1
      do 3 np = 6 , 9
        if(np.eq.7) go to 3
       do 2 k = 1 , nh
        snp=0.
        ssp=0.
        do 1 j = 1 , ids
         snp=snp+pgl(np,k,2,j)
         ssp=ssp+pgl(np,k,i2,j)
  1     continue
        unp=snp/ids
        usp=ssp/ids
        do 5 j = 1 , ids
          pgl(np,k,1,j)=unp
          pgl(np,k,its,j)=usp
  5     continue
  2    continue
  3   continue
      return
      end

