!  ver jan2020 2-D eddy diffusion coefficient
!  quasi-uniform grid
!  O2 в приближении малой компоненты  (прогонка)
      subroutine o2_eddy (an11,an1,an2,an3,an6,vr,vi,vj,
     *                  q,eddyco,rads,rp,g,n1,n2,n,dt)

      USE mo_bas_gsm, ONLY:amo2,amn2,amo,bk,pi,re 

      dimension an1(n1,n2,n),an2(n1,n2,n),
     *          an3(n1,n2,n),an6(n1,n2,n),an11(n1,n2,n),
     *          vr(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n),
     *          q(n1,n2,n),rads(n),rp(n),g(n),eddyco(n,n1)
      allocatable a(:),b(:),c(:),f(:),cmd(:)
     *         ,h(:),alf(:),bet(:),hsr(:),co2(:)

      data gam/1.e-20/,geom/1.1/

      allocate (a(n),b(n),c(n),f(n),cmd(n)
     *         ,h(n),alf(n),bet(n),hsr(n),cO2(n))
      const=bk/amO2
c*******
      dtet=pi/(n1-1)
      dfi=2.*pi/n2
      do  i=2,n1-1
        teta=dtet*(i-1)
        sin_t=sin(teta)
        cot_t=cos(teta)/sin_t
        do j=1,n2
         jp=j+1
         jm=j-1
         if(j.eq.n2) jp=1
         if(j.eq.1)  jm=n2
         cO2(1)=an1(i,j,1)
         do k=1,n
           h(k)=const*an6(i,j,k)/g(k)
           sum=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
           ams=(amO2*an1(i,j,k)+amN2*an2(i,j,k)+amO*an3(i,j,k))/sum
           hsr(k)=bk*an6(i,j,k)/(ams*g(k))
c    . . . Coef. Mol. Dif.
c         epok=1
           epok=exp(2.8/an6(i,j,k))
           sum1=an1(i,j,k)+an2(i,j,k)
           sum2=an1(i,j,k)+an3(i,j,k)
	   if(sum1.le.0.) then

	print*,'o2pro',i,j,k, an1(i,j,k-1),an2(i,j,k-1),an3(i,j,k-1)
!	pause
             an1(i,j,k)=an1(i,j,k-1)
	     an2(i,j,k)= an2(i,j,k-1)
	     an3(i,j,k)=an3(i,j,k-1)
	     sum1=an1(i,j,k)+an2(i,j,k)
             sum2=an1(i,j,k)+an3(i,j,k)

	   end if
           d12=0.829e17/sum1*an6(i,j,k)**0.724*epok
           d13=0.969e17/sum2*an6(i,j,k)**0.774*epok
          obr=(an2(i,j,k)/d12+an3(i,j,k)/d13)/(an2(i,j,k)+an3(i,j,k))
          cmd(k)=1./obr
! test
!          cmd(k)=3.e17/sum*sqrt(an6(i,j,k))  ! к-т м.диффузии NO
! test end
          alf(k)=cmd(k)/(cmd(k)+eddyco(k,i))
          bet(k)=eddyco(k,i)/(cmd(k)+eddyco(k,i))
         end do
         c o2(1)=an1(i,j,1)
!  for uniform grid
         con1=alog(geom)/(geom-1.) 
         do k=2,n-1
          rk=rads(k)+re
          pro= (rp(k)+rp(k-1))*.5
          c o2(k)=an1(i,j,k)
          sum=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
c    . . . источники и потери
          alpha=9.9e-34*exp(470./an6(i,j,k))
          betta=1.1e-34*exp(510./an6(i,j,k))
          ! if(q(i,j,k).ne.0.) betta=0.
c    . . .
          qq=(alpha*sum+gam)*an3(i,j,k)**2
          p=q(i,j,k)+betta*an3(i,j,k)*sum
c ============================================
c          p=0.
c          qq=0.
          dtnp=alog(an6(i,j,k+1)/an6(i,j,k))
          dtnm=alog(an6(i,j,k)/an6(i,j,k-1))
c      . . .  к-т диффузии в дробной точке с учетом неравном сетки
          cmdp=(cmd(k+1)+cmd(k)+eddyco(k+1,i)+eddyco(k,i))*.5
          cmdm=(cmd(k)+cmd(k-1)+eddyco(k,i)+eddyco(k-1,i))*.5
        
          cMolp=(cmd(k+1)+geom*cmd(k))/(geom+1.)
          cMolm=(cmd(k)+geom*cmd(k-1))/(geom+1.)
          cTdp=(eddyco(k+1,i)+geom*eddyco(k,i))/(geom+1.)
          cTdm=(eddyco(k,i)+geom*eddyco(k-1,i))/(geom+1.)
          cAllp=cMolp+cTdp
          cAllm=cMolm+cTdm

          a(k)=(cAllp/pro)/(rp(k)*con1)
          c(k)=(cAllm/pro)/(rp(k-1)*con1)
          b(k)=a(k)+c(k)+1./dt+p
          hkp=(h(k+1)+geom*h(k))/(geom+1.)
          hkm=(h(k)+geom*h(k-1))/(geom+1.)
          adop=cMolp*0.5*(1./hkp+dtnp/(rp(k)*con1))/pro
          cdop=-cMolm*0.5*(1./hkm+dtnm/(rp(k-1)*con1))/pro
          a(k)=a(k)+adop
          c(k)=c(k)+cdop
          b(k)=b(k)-adop-cdop
          hsrp=(hsr(k+1)+geom*hsr(k))/(geom+1.)
          hsrm=(hsr(k)+geom*hsr(k-1))/(geom+1.)
          aturb=cTdp*0.5*(1./hsrp+dtnp/(rp(k)*con1))/pro
          cturb=-cTdm*0.5*(1./hsrm+dtnm/(rp(k-1)*con1))/pro
          a(k)=a(k)+aturb
          c(k)=c(k)+cturb
          b(k)=b(k)-aturb-cturb
! d(nv)/dr -> vp*(n(k+1)-n(k)) + vm(n(k)-n(k-1))
          vrp=0.5*(vr(i,j,k)+abs(vr(i,j,k)))
          vrm=0.5*(vr(i,j,k)-abs(vr(i,j,k)))
          a(k)=a(k)-vrm/(rp(k)*con1)
          c(k)=c(k)+vrp/(rp(k-1)*con1)
          b(k)=b(k)-vrm/(rp(k)*con1)+vrp/(rp(k-1)*con1)
          b(k)=b(k)+0.5*(1.+sign(1.,vr(i,j,k)))*
     *              (vr(i,j,k+1)-vr(i,j,k))/(rp(k)*con1)
     *             +0.5*(1.-sign(1.,vr(i,j,k)))*
     *              (vr(i,j,k)-vr(i,j,k-1))/(rp(k-1)*con1)

!!!!!!!!!!!!! end
          f(k)=qq+co2(k)/dt
c . . .  Дивергенция V:
          del=2.*dfi*rk*sin_t
          div=(vj(i,jp,k)-vj(i,jm,k))/del
c . . .  учет котангенса:
          del=2.*dtet*rk
          div=div+vi(i,j,k)/rk*cot_t+
     *         (vi(i+1,j,k)-vi(i-1,j,k))/del
          if(div.gt.0.) then
             b(k)=b(k)+div
          else
             f(k)=f(k)-div*co2(k)
          end if
        end do
c . . . элемент массива используется для удобства
c        f(n)=h(n)/(rp(n-1)+h(n))
c        f(n)=exp(-rp(n-1)/h(n))
        obrh=(rp(n)/h(n)+geom*rp(n-1)/h(n-1))*con1/(geom+1.)
        f(n)=exp(-obrh)
        kiss=2
        call progon (co2,a,b,c,f,n,kiss)
c     . . .
! . . . reASSIGN 

        an11(i,j,1)=co2(1)
!!        do k=2,n
!!          an11(i,j,k)=co2(k)
!!        end do
!!	...  and correction
        do k=n-1,2,-1
         if(co2(k+1).gt.co2(k)) then 
	 ! profil correction 31/03/15
           co2(k)=co2(k+1)*exp(0.5*(rp(k)/h(k+1)+rp(k-1)/h(k)))
	   print*, 'nO2 correction!!!', co2(k),i,j,k  
         end if	  
         an11(i,j,k)=c o2(k)
        end do
       end do
      end do
	deallocate (a,b,c,f,cmd
     *         ,h,alf,bet,hsr,c o2)
      return
      end
c    . . . Прогонка
      subroutine progon(dim,a,b,c,f,n,kiss)
      dimension dim(n),a(n),b(n),c(n),f(n)
     *         ,pa(100),pb(100)
      pa(2)=0.
      pb(2)=dim(1)
      n1=n-1
      do 1 k=2,n1
        del=b(k)-pa(k)*c(k)
        pa(k+1)=a(k)/del
        pb(k+1)=(c(k)*pb(k)+f(k))/del
c        if(pa(k+1).gt.1) then
c         print 100,pa(k+1),k
c100     format(' pa.gt.1',e8.1,i4)
c      end if
   1  continue
      if(kiss.eq.1) dim(n)=pb(n)/(1.-pa(n))           ! for temperature
      if(kiss.eq.2) dim(n)=f(n)*pb(n)/(1.-pa(n)*f(n)) ! for composition
      do 2 l=2,n
        k=n-l+1
        dim(k)=pa(k+1)*dim(k+1)+pb(k+1)
    2 continue
      return
      end
