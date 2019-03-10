      subroutine bongi(gins,ins,nh,its,ids)
      dimension gins(ins,nh,its,ids)
      i2=its-1
      do 3 i=1,3
          np=i
          if(i.eq.2)np=5
          if(i.eq.3)np=6
       do 2 k = 1 , nh
c      ssp - sum s.pole
c      snp - sum n.pole
        s np=0.
        s sp=0.
        do 1 j = 1 , ids
         snp=snp+gins(np,k,2,j)
         ssp=ssp+gins(np,k,i2,j)
   1    continue
c
        unp=snp/ids
        usp=ssp/ids
        do 4 j=1,ids
          gins(np,k,1,j)=unp
          gins(np,k,its,j)=usp
    4   continue
    2  continue
    3 continue
      return
      end
