       subroutine connot(pgl,rads,kpars,nh,its,ids)
c      . . . расчет NO по Куликову
       dimension pgl(kpars,nh,its,ids),rads(nh)
       do 1 k = 1 , nh
        z=rads(k)*1.e-5
        pok=z-110.
        pok=pok*pok/900.
        cno=exp(-pok)*2.7e7
        pok1=(z-140.)/3.29*100.
        do 2 j = 1 , ids
         do 3 i = 1 , its
          if(k.le.12) then
           pgl(4,k,i,j)=cno
          else
c          pgl(4,k,i,j)=1.e7*exp(-pok1/pgl(7,nh,i,j))
           pgl(4,k,i,j)=2.e7*exp(-pok1/pgl(7,nh,i,j))
          endif
    3    continue
    2   continue
    1  continue
       print *,' connot - end'
       return
       end
