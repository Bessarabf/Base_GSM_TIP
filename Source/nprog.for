      subroutine nprog(cNd,cNoi,cO2i,cNe
     *                ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
c    N(1s) progonka as O2pro
      dimension pgl(kpars,nh,its,ids),cNd(its,ids,nh),
     *          cO2i(nh),cNOi(nh),cNe(nh),
     *          rads(nh),rp(nh),g(nh),ctd(nh)
    
      allocatable a(:),b(:),c(:),f(:),cmd(:),cN(:),
     *            h(:),hsr(:),alf(:),bet(:)
      INCLUDE 'alpha.inc'	 
      data amo2,amn2,amo,amn/ 53.12e-24,46.51e-24,26.56e-24,23.26e-24/
     *    ,bk/1.38e-16/,re/6.371e8/,pi/3.14159/
      allocate (a(nh),b(nh),c(nh),f(nh),cmd(nh),cN(nh),
     *            h(nh),hsr(nh),alf(nh),bet(nh))

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
      cN(1)=pgl(5,1,i,j)
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
          cmd(k)=4.55e17/sum*sqrt(pgl(7,k,i,j))  ! к-т м.диффузии N
          alf(k)=cmd(k)/(cmd(k)+ctd(k))
          bet(k)=ctd(k)/(cmd(k)+ctd(k))
         end do
         do k=2,nh-1
          rk=rads(k)+re
          pro= (rp(k)+rp(k-1))*.5
          cN(k)=pgl(5,k,i,j)
          sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
          ot=300./pgl(9,k,i,j)
          ots=sqrt(ot)

c p -lost, and q - source
          u=alfa1*ot**0.85*(1.-r1)*cNoi(k)*cNe(k)
          u=u+(alyam9*pgl(3,k,i,j)+alyam11/ots*cNe(k)+
     *    alyam16*pgl(2,k,i,j)+alyam13)*cNd(i,j,k)
          q=u+r3*alyam6*pgl(14,k,i,j)
          ho2=const1*pgl(7,k,i,j)/g(k)
          x=(rads(k)+re)/h o2
          ch=chept(x,hi)
          pok=1.e-8*(pgl(1,k,i,j)*h o2*ch*1.E-4)**0.38
          if (pok.gt.30.) pok=30.
          q=q+alyam14*exp(-pok)*pgl(4,k,i,j)
c . . .  lost
          w=alyam2*cO2i(k)
          w=w+alyam8*exp(-3220./pgl(7,k,i,j))*pgl(1,k,i,j)
          p=w+alyam10*pgl(4,k,i,j)
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
          a(k)=a(k)+aprim-pgl(10,k+1,i,j)/pro*.5
          c(k)=c(k)+cprim+pgl(10,k-1,i,j)/pro*.5
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
         pgl(5,k,i,j)=cN(k)
        end do
      deallocate (a ,b ,c ,f ,cmd ,cN ,
     *            h ,hsr ,alf ,bet )

      return
      end
