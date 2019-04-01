      subroutine noprog(cNd,cNo,cO2i
     *                 ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
c    NO в приближении малой компоненты  (прогонка)
      dimension pgl(kpars,nh,its,ids),
     *          cO2i(nh),cNo(nh),cNd(its,ids,nh),
     *          rads(nh),rp(nh),g(nh),ctd(nh)
      allocatable a(:),b(:),c(:),f(:),cmd(:)
     *          ,h(:),hsr(:),alf(:),bet(:)
	allocate (a(nh),b(nh),c(nh),f(nh),cmd(nh)
     *          ,h(nh),hsr(nh),alf(nh),bet(nh))

      INCLUDE 'alpha.inc'	 
      data amo2,amn2,amo,amno/ 53.12e-24,46.51e-24,26.56e-24,49.82e-24/
     *    ,bk/1.38e-16/,re/6.371e8/,pi/3.14159/
      const=bk/amno
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
      cNO(1)=pgl(4,1,i,j)
      do k=1,nh     
c     . . . scale height
         h(k) =const*pgl(7,k,i,j)/g(k)
         sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
         ams=(amo*pgl(3,k,i,j)+amo2*pgl(1,k,i,j)+
     *        amn2*pgl(2,k,i,j))/sum
         hsr(k)=bk*pgl(7,k,i,j)/(ams*g(k))
c    . . . Coef. Mol. Dif.
c         epok=1
!          epok=exp(2.8/an6(i,j,k))
!          sum1=an1(i,j,k)+an2(i,j,k)
!          sum2=an1(i,j,k)+an3(i,j,k)
!          d12=0.829e17/sum1*an6(i,j,k)**0.724*epok
!          d13=0.969e17/sum2*an6(i,j,k)**0.774*epok
!          obr=(an2(i,j,k)/d12+an3(i,j,k)/d13)/(an2(i,j,k)+an3(i,j,k))
!          cmd(k)=1./obr
          cmd(k)=3.e17/sum*sqrt(pgl(7,k,i,j))  ! к-т м.диффузии NO
          alf(k)=cmd(k)/(cmd(k)+ctd(k))
          bet(k)=ctd(k)/(cmd(k)+ctd(k))
         end do
         do k=2,nh-1
          rk=rads(k)+re
          pro= (rp(k)+rp(k-1))*.5
          cNO(k)=pgl(4,k,i,j)
          sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
! . .  q - sourse and...
         u=(alyam7*cNd(i,j,k)+
     *      alyam8*exp(-3220./pgl(7,k,i,j))*pgl(5,k,i,j))*pgl(1,k,i,j)
         q=u+alyam3*cO2i(k)*pgl(2,k,i,j)
! . . . p - lost
         w=alyam10*pgl(5,k,i,j)+alyam12*cNd(i,j,k)
         w=w+alyam15*pgl(3,k,i,j)*pgl(2,k,i,j)*exp(940./pgl(7,k,i,j))
 !       p=w+alyam1*cO2i(k)
         p=w+alyam1*cO2i(k)+pgl(15,k,i,j)/pgl(4,k,i,j)
 	 
         ho2=const1*pgl(7,k,i,j)/g(k)
         x=(rads(k)+re)/h o2
         ch=chept(x,hi)
         pok=1.e-8*(pgl(1,k,i,j)*h o2*ch*1.e-4)**0.38
         if (pok.gt.30.) pok=30.
         p=p+alyam14*exp(-pok)
c ============================================
c          p=0.
c          qq=0.
          dtnp=alog(pgl(7,k+1,i,j)/pgl(7,k,i,j))
          dtnm=alog(pgl(7,k,i,j)/pgl(7,k-1,i,j))
c      . . . к-т диффузии в дробной точке
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
!          a(k)=a(k)+aprim-pgl(10,k+1,i,j)/pro*.5
!          c(k)=c(k)+cprim+pgl(10,k-1,i,j)/pro*.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!30.10.18 !!!!!!!!!!!!!!!!!!!!
          a(k)=a(k)+aprim !30.10.18
          c(k)=c(k)+cprim !30.10.18
          if(pgl(10,k,i,j).le.0.) then
           a(k)=a(k)-pgl(10,k+1,i,j)/rp(k)
           b(k)=b(k)-pgl(10,k,i,j)/rp(k)
          else  
           b(k)=b(k)+pgl(10,k,i,j)/rp(k-1)
           c(k)=c(k)+pgl(10,k-1,i,j)/rp(k-1)
          end if 
!!!!!!!!!!!!!!!!!!!!! end 30.10.18 !!!!!!!!!!!!!!!!!!!!!
          b(k)=b(k)+cl*(cmd(k)+ctd(k))/pro
          f(k)=q+cNO(k)/dt
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
             f(k)=f(k)-div*cNO(k)
          end if
        end do
c . . . элемент массива используется для удобства
c       
        f(nh)=exp(-rp(nh-1)/h(nh))
        kiss=2
        call progon (cNO,a,b,c,f,nh,kiss)
c     . . .
        do k=1,nh
         pgl(4,k,i,j)=cNO(k)
        end do
      deallocate (a,b,c,f,cmd
     *          ,h,hsr,alf,bet)
      return
      end
