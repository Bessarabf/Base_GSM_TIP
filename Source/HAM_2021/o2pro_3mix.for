!
!  3-composition mixture                        23.12.2020
!  ver jan2020 2-D eddy diffusion coefficient   
!  NON - quasi-uniform grid                     21.12.2020 
!  
      subroutine o2pro_3mix(an11,an1,an2,an3,an6,vr,vi,vj,
     *                  q,eddyco,ro,rads,rp,g,n1,n2,n,dt)

      USE mo_gsm_const, ONLY:amo2,amn2,amo,bk,pi,re

      dimension an1(n1,n2,n),an2(n1,n2,n),
     *          an3(n1,n2,n),an6(n1,n2,n),an11(n1,n2,n),
     *          vr(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n),ro(n1,n2,n),
     *          q(n1,n2,n),rads(n),rp(n),g(n),eddyco(n,n1)
     *          ,dp(3,3)
      allocatable a(:),b(:),c(:),f(:),cmd(:)
     *         ,h(:),alf(:),bet(:),hsr(:),cO2(:)
     *         ,dp12(:),dp13(:) 
      
      data gam/1.e-20/

      allocate (a(n),b(n),c(n),f(n),cmd(n)
     *         ,h(n),alf(n),bet(n),hsr(n),cO2(n),
     *          dp12(n),dp13(n))
      const=bk/amO2
      constN2=bk/amN2
      constO =bk/amO
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
           cono2=an1(i,j,k)
           conn2=an2(i,j,k)
           cono=an3(i,j,k)
           tempp=an6(i,j,k)
           rog=ro(i,j,k)
           call kmoldn(dp,cono2,conn2,cono,tempp,rog)
           cmd(k)=dp(1,1)
           ! diffusion O2 - N2
           dp12(k)=dp(1,2)
           ! diffusion O2 - O
           dp13(k)=dp(1,3)
!!! test 
!!         cmd(k)=3.e17/sum*sqrt(an6(i,j,k))  ! k-t mol dif for NO
!!!
           alf(k)=cmd(k)/(cmd(k)+eddyco(k,i))
           bet(k)=eddyco(k,i)/(cmd(k)+eddyco(k,i))
         end do
         cO2(1)=an1(i,j,1)
         do k=2,n-1
          rk=rads(k)+re
          pro= (rp(k)+rp(k-1))*.5
          c o2(k)=an1(i,j,k)
          sum=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
c    . . . источники и потери
          alpha=9.9e-34*exp(470./an6(i,j,k))
          betta=1.1e-34*exp(510./an6(i,j,k))
          if(q(i,j,k).ne.0.)betta=0.
c    . . .
          qq=(alpha*sum+gam)*an3(i,j,k)**2
          p=q(i,j,k)+betta*an3(i,j,k)*sum
c ============================================
c          p=0.
c          qq=0.
          dtnp=alog(an6(i,j,k+1)/an6(i,j,k))
          dtnm=alog(an6(i,j,k)/an6(i,j,k-1))
c      . . . к-т диффузии в дробной точке
          cmdp=(cmd(k+1)+cmd(k)+eddyco(k+1,i)+eddyco(k,i))*.5
          cmdm=(cmd(k)+cmd(k-1)+eddyco(k,i)+eddyco(k-1,i))*.5
c          vrp=(vr(i,j,k)+abs(vr(i,j,k)))*.5
c          vrm=(vr(i,j,k)-abs(vr(i,j,k)))*.5
c          a(k)=(cmdp/pro-vrm)/rp(k)
c          c(k)=(cmdm/pro+vrp)/rp(k-1)
          a(k)=(cmdp/pro)/rp(k)
          c(k)=(cmdm/pro)/rp(k-1)
          b(k)=a(k)+c(k)+1./dt+p
          clp=0.5*(alf(k+1)/h(k+1)+dtnp/rp(k)+bet(k+1)/hsr(k+1))
          clm=-0.5*(alf(k-1)/h(k-1)+dtnm/rp(k-1)+bet(k-1)/hsr(k-1))
          aprim=(cmd(k+1)+eddyco(k+1,i))/pro*clp
          cprim=(cmd(k-1)+eddyco(k-1,i))/pro*clm
          cl=0.5*(-dtnp/rp(k)+dtnm/rp(k-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!30.10.18 !!!!!!!!!!!!!!!!!!!!
          a(k)=a(k)+aprim !30.10.18
          c(k)=c(k)+cprim !30.10.18
          if(vr(i,j,k).le.0.) then
           a(k)=a(k)-vr(i,j,k+1)/rp(k)
           b(k)=b(k)-vr(i,j,k)/rp(k)
          else  
           b(k)=b(k)+vr(i,j,k)/rp(k-1)
           c(k)=c(k)+vr(i,j,k-1)/rp(k-1)
          end if 
!!!!!!!!!!!!!!!!!!!!! end 30.10.18 !!!!!!!!!!!!!!!!!!!!!

          b(k)=b(k)+cl*(cmd(k)+eddyco(k,i))/pro
          f(k)=qq+co2(k)/dt
!!!!!!!!!!  23.12.2020 O2 - N2 diffusion
          hN2p=constN2*an6(i,j,k+1)/g(k+1)
          hN2=constN2*an6(i,j,k)/g(k)
          hN2m=constN2*an6(i,j,k-1)/g(k-1)
          aN2p=an2(i,j,k+1)
          aN2m=an2(i,j,k-1)
          difP=(an2p-an2(i,j,k))/rp(k)
          difM=(an2(i,j,k)-aN2m)/rp(k-1)
       ! k+1/2 point 
          fluxP=difP+0.5*(aN2p+an2(i,j,k))*(0.5*(hN2p+hN2)/(hN2p*hN2)+
     *          dtnp/rp(k))
          fluxP=fluxP*(dp12(k+1)+dp12(k))
       ! k-1/2 point
          fluxM=difM+0.5*(aN2m+an2(i,j,k))*(0.5*(hN2m+hN2)/(hN2m*hN2)+
     *          dtnm/rp(k-1)) 
          fluxM=fluxM*(dp12(k)+dp12(k-1))
        !  add to f(k)
          f(k)=f(k)+(fluxP-fluxM)/pro
 
!!!!!!!!!!  23.12.2020 O2 - O diffusion
          hOp=constO*an6(i,j,k+1)/g(k+1)
          hO=constO*an6(i,j,k)/g(k)
          hOm=constO*an6(i,j,k-1)/g(k-1)
          aOp=an3(i,j,k+1)
          aOm=an3(i,j,k-1)
          difPo=(aOp-an3(i,j,k))/rp(k)
          difMo=(an3(i,j,k)-aOm)/rp(k-1)
       ! k+1/2 point 
          fluxPo=difPo+0.5*(aOp+an3(i,j,k))*(0.5*(hOp+hO)/(hOp*hO)+
     *          dtnp/rp(k))
          fluxPo=fluxPo*(dp13(k+1)+dp13(k))
       ! k-1/2 point
          fluxMo=difMo+0.5*(aOm+an3(i,j,k))*(0.5*(hOm+hO)/(hOm*hO)+
     *          dtnm/rp(k-1)) 
          fluxMo=fluxMo*(dp13(k)+dp13(k-1))
        !  add to f(k)
          f(k)=f(k)+(fluxPo-fluxMo)/pro      
         
c . . .  div V :
          del=2.*dfi*rk*sin_t
          div=(vj(i,jp,k)-vj(i,jm,k))/del
c . . .  cotan :
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
        f(n)=exp(-rp(n-1)/h(n))
        kiss=2
        call progon (c o2,a,b,c,f,n,kiss)
c     . . .
! . . . reASSIGN 

        an11(i,j,1)=cO2(1)
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
     *         ,h,alf,bet,hsr,cO2,dp12,dp13)
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
