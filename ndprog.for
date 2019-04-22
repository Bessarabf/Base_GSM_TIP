      subroutine ndnew(cNd,cNo,cNoi,cN2i,cNe
     *                ,pgl,kpars,nh,its,ids,i,j,key,dt)
      dimension cNd(its,ids,nh),cNoi(nh),cN2i(nh),cNe(nh),
     *          cNo(nh),pgl(kpars,nh,its,ids)
      INCLUDE 'alpha.inc'
c
c     data  alfa1/4.2e-7/,r1,r3/0.7,0.5/
c    *    , alyam4 /1.4e-10/
c    *    , alyam6,  alyam7,   alyam9
c    *    /    0.1,  5.e-12,  4.5e-13/
c    *    ,alyam11, alyam12,  alyam13, alyam16
c    *    /3.6e-10,  7.e-11,  1.06e-5, 2.3e-14/

      do 1 k=1,nh
       ot=300./pgl(9,k,i,j)
       ots=sqrt(ot)
       tr=(pgl(8,k,i,j)+pgl(7,k,i,j))*.5
c p -lost, and q- source
c
       a=alfa1*ot**0.85*r1*cNoi(k)*cNe(k)
       a=a+alyam4*(300./tr)**0.44*cN2i(k)*pgl(3,k,i,j)
       q=a+r3*alyam6*pgl(14,k,i,j)
       b=alyam7*pgl(1,k,i,j)
       b=b+alyam9*pgl(3,k,i,j)
       b=b+alyam11/ots*cNe(k)
       b=b+alyam12*cNo(k)+alyam16*pgl(2,k,i,j)
       p=b+alyam13
      ! if(key.eq.0) then
      !  cNd(i,j,k)=q/p
      ! else
        cNd(i,j,k)=(cNd(i,j,k)+q*dt)/(1.+p*dt)
      ! end if
   1  continue
      return
      end
