      subroutine ndprog(cNd,cNo,cNoi,cN2i,cNe
     *                ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
       dimension pgl(kpars,nh,its,ids),cNd(its,ids,nh),
     *          cNO(nh),cNOi(nh),cNe(nh),cN2i(nh),
     *          rads(nh),rp(nh),g(nh),ctd(nh)

c    N(2D) progonka as O2pro
      allocatable a(:),b(:),c(:),f(:),cmd(:),cN(:),
     *            h(:),hsr(:),alf(:),bet(:)


      INCLUDE 'alpha.inc'
      data amo2,amn2,amo,amn/ 53.12e-24,46.51e-24,26.56e-24,23.26e-24/
     *    ,bk/1.38e-16/,re/6.371e8/,pi/3.14159/
      allocate (a(nh),b(nh),c(nh),f(nh),cmd(nh),cN(nh),
     *            h(nh),hsr(nh),alf(nh),bet(nh))

c
c     data  alfa1/4.2e-7/,r1,r3/0.7,0.5/
c    *    , alyam4 /1.4e-10/
c    *    , alyam6,  alyam7,   alyam9
c    *    /    0.1,  5.e-12,  4.5e-13/
c    *    ,alyam11, alyam12,  alyam13, alyam16
c    *    /3.6e-10,  7.e-11,  1.06e-5, 2.3e-14/
      const=bk/amn
      const1=bk/amo2
c*******
      dtet=pi/(its-1)
      dfi=2.*pi/ids
      teta=dtet*(i-1)
      sin_t=sin(teta)
      cot_t=cos(teta)/sin_t
      jp=j+1
      jm=j-1
      if(j.eq.ids) jp=1
      if(j.eq.1)  jm=ids
      cN(1)=cNd(i,j,1)
      do k=1,nh     
c     . . . scale heights
         h(k) =const*pgl(7,k,i,j)/g(k)
         sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
         ams=(amo*pgl(3,k,i,j)+amo2*pgl(1,k,i,j)+
     *        amn2*pgl(2,k,i,j))/sum
         hsr(k)=bk*pgl(7,k,i,j)/(ams*g(k))
c     . . .
         cmd(k)=4.55e17/sum*sqrt(pgl(7,k,i,j))  ! mol dif N
         alf(k)=cmd(k)/(cmd(k)+ctd(k))
         bet(k)=ctd(k)/(cmd(k)+ctd(k))
      end do
      do k=2,nh-1
         ot=300./pgl(9,k,i,j)
         ots=sqrt(ot)
         tr=(pgl(8,k,i,j)+pgl(7,k,i,j))*.5
         rk=rads(k)+re
         pro=(rp(k)+rp(k-1))*.5
         cN(k)=cNd(i,j,k)
         sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
         ot=300./pgl(9,k,i,j)
         ots=sqrt(ot)
c . . . p -lost, and q- source
         q=alfa1*ot**0.85*r1*cNoi(k)*cNe(k)               ! NO+ + e
         q=q+alfa3*ot**0.4*r2*cN2i(k)*cNe(k)              ! N2+ + e
         q=q+alyam4*(300./tr)**0.44*cN2i(k)*pgl(3,k,i,j)  ! N2+ + O
         q=q+r3*alyam6*pgl(14,k,i,j)                      ! q(N2+) - dissosiation

         p=alyam7*pgl(1,k,i,j)
         p=p+alyam9*pgl(3,k,i,j)
         p=p+alyam11/ots*cNe(k)
         p=p+alyam12*cNo(k)+alyam16*pgl(2,k,i,j)
         p=p+alyam13
c . . .
         dtnp=alog(pgl(7,k+1,i,j)/pgl(7,k,i,j))
         dtnm=alog(pgl(7,k,i,j)/pgl(7,k-1,i,j))
c      . . . dif in semi integer point
         cmdp=(cmd(k+1)+cmd(k)+ctd(k+1)+ctd(k))*.5
         cmdm=(cmd(k)+cmd(k-1)+ctd(k)+ctd(k-1))*.5

         a(k)=(cmdp/pro)/rp(k)
         c(k)=(cmdm/pro)/rp(k-1)
         b(k)=a(k)+c(k)+1./dt+p
         clp=0.5*(alf(k+1)/h(k+1)+dtnp/rp(k)+bet(k+1)/hsr(k+1))
         clm=-0.5*(alf(k-1)/h(k-1)+dtnm/rp(k-1)+bet(k-1)/hsr(k-1))
         aprim=(cmd(k+1)+ctd(k+1))/pro*clp
         cprim=(cmd(k-1)+ctd(k-1))/pro*clm
         cl=0.5*(-dtnp/rp(k)+dtnm/rp(k-1))
!          a(k)=a(k)+aprim-pgl(10,k+1,i,j)/pro*.5 ! var before 30.10.18
!          c(k)=c(k)+cprim+pgl(10,k-1,i,j)/pro*.5 !
!!! like O2pro 30/10/18 !!!!!!!!!!!!!!!!!!!!!!
          a(k)=a(k)+aprim
          c(k)=c(k)+cprim
          if(pgl(10,k+1,i,j).le.0.) then
           a(k)=a(k)-pgl(10,k+1,i,j)/rp(k)
           b(k)=b(k)-pgl(10,k,i,j)/rp(k)
          else  
           b(k)=b(k)+pgl(10,k,i,j)/rp(k-1)
           c(k)=c(k)+pgl(10,k-1,i,j)/rp(k-1)
          end if 
!!!!!!!!!!!!!!!!!!!!! end 30.10.18 !!!!!!!!!!!!!!!!!!!!!
          b(k)=b(k)+cl*(cmd(k)+ctd(k))/pro
          f(k)=q+cN(k)/dt
c . . .  Дивергенция V:
          del=2.*dfi*rk*sin_t
          div=(pgl(12,k,i,jp)-pgl(12,k,i,jm))/del
c . . .  учет котангенса:
          del=2.*dtet*rk
          div=div+pgl(11,k,i,j)/rk*cot_t+
     *         (pgl(11,k,i+1,j)-pgl(11,k,i-1,j))/del
          if(div.gt.0.) then
             b(k)=b(k)+div
          else
             f(k)=f(k)-div*cN(k)
          end if
        end do
c . . . элемент массива используется для удобства
c       
        f(nh)=exp(-rp(nh-1)/h(nh))
        kiss=2
        call progon (cN,a,b,c,f,nh,kiss)
c     . . .
      do k=1,nh
         cNd(i,j,k)=cN(k)
      end do
      deallocate (a,b,c,f,cmd,cN,h,hsr,alf,bet )
      return
      end
