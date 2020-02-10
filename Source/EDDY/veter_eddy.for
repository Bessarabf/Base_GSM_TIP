c veter
c vi_r
c vj_r
c vn_t
c vn_fi
c vrprim_s
! 
! version jan.2020 with 2-D eddy diffusion coefficient

      subroutine veter_eddy(vi1,vj1,vi,vj,vr,vim,vjm,an1,an2,an3,an6,ro,
     *                 pgl,pgi,rp,r,eddyco,n,n1,n2,kpars,ins,dt)
!

      dimension vi(n1,n2,n),vi1(n1,n2,n),vj(n1,n2,n),vj1(n1,n2,n),
     *          vr(n1,n2,n),vim(n,n1,n2),vjm(n,n1,n2),
     *          an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n),
     *          an6(n1,n2,n),ro(n1,n2,n),
     *          pgl(kpars,n,n1,n2),pgi(ins,n,n1,n2),
     *          rp(n),r(n),eddyco(n,n1)
c     . . .   Расчет меридиональной компоненты
c     . . .  Следующий порядок расщепления не произволен!!!
c     . . . Расщепление вдоль меридиана
      key=0 ! менять знак для скор. при переходе через полюс
      call vn_t(vi1,vi,r,n,n1,n2,dt,key)
c     . . . Расщепление вдоль азимута
      call vn_fi(vi1,vj,r,n,n1,n2,dt)
c     . . . Расщепление по вертикали
      call vi_r_ham(vi1,vi,vj,vr,vim,an1,an2,an3,an6,ro,
     *          pgl,pgi,rp,r,eddyco,n,n1,n2,kpars,ins,dt)
c     . . .   Расчет азимутальной компоненты
c     . . .  Следующий порядок расщепления не произволен!!!
c     . . . Расщепление вдоль меридиана
!      key=0
      key=1 ! не меняется знак скорости при переходе через полюс
      call vn_t(vj1,vi,r,n,n1,n2,dt,key)
c     . . . Расщепление вдоль азимута
      call vn_fi(vj1,vj,r,n,n1,n2,dt)
c     . . . Расщепление по вертикали
      call vj_r_ham(vj1,vj,vi,vr,vjm,an1,an2,an3,an6,ro,
     *          pgl,pgi,rp,r,eddyco,n,n1,n2,kpars,ins,dt)
      return
      end
      subroutine vi_r_ham(vi1,vi,vj,vr,vim,an1,an2,an3,an6,ro,
     *                pgl,pgi,rp,r,eddyco,n,n1,n2,kpars,ins,dt)
      dimension vi(n1,n2,n),vi1(n1,n2,n),vj(n1,n2,n),
     *          vr(n1,n2,n),vim(n,n1,n2),
     *          an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n),
     *          an6(n1,n2,n),ro(n1,n2,n),
     *          pgl(kpars,n,n1,n2),pgi(ins,n,n1,n2),
     *          dragViGSM(n1,n2,n),
     *          rp(n),r(n),eddyco(n,n1),
     *          pa(31),pb(31),amu(31)
      data am1,am2,am3,bk,om,pi/53.12e-24,46.51e-24,26.56e-24,
     *     1.38e-16,7.27e-5,3.14159/,re/6.371e08/
      data rg/8.31442e07/
c     . . . средняя масса
      am_s(a1,a2,a3)=(32.*a1+28.*a2+16.*a3)/
     *               (a1+a2+a3)
c     . . . вычисление шагов по долготе и широте
      dfi=2.*pi/n2
      n1m1=n1-1
      dtet=pi/(n1m1)
      pa(2)=0.
      do i=2,n1m1
c
        tet=dtet*(i-1)
        cos_t=cos(tet)
        sin_t=sin(tet)
        do j=1,n2
          costet=.98*cos_t-.196*sin_t*cos(dfi*(j-1))
          pb(2)=vi(i,j,1)
          do k=1,n
            amu(k)=3.34e-6*(an6(i,j,k)**0.71)
            cp=3.5*(an1(i,j,k)+an2(i,j,k))+2.5*an3(i,j,k)
            cp=cp*bk/(am1*an1(i,j,k)+am2*an2(i,j,k)+
     *         am3*an3(i,j,k))
c       . . . Турбулентная вязкость
            amut=2.*eddyco(k,i)/cp
            amu(k)=amu(k)+amut
          end do
          do k=2,n-1
            rk=r(k)+re
            ot=dt/rp(k)
            del=2.*rp(k)*ro(i,j,k)
            a=dt*(amu(k+1)/rp(k+1)+amu(k)/rp(k))/del
            c=dt*(amu(k)/rp(k)+amu(k-1)/rp(k-1))/del
            vp=0.5*(vr(i,j,k)+abs(vr(i,j,k)))
            vm=0.5*(vr(i,j,k)-abs(vr(i,j,k)))
            vm=vm*ot
            vp=vp*ot
            a=a-vm
            c=c+vp
            b=1.+a+c
c . . .    правая часть
c            f=vi(i,j,k)
c . . .   для выбранного варианта расщепления
            f=vi1(i,j,k)
            f=f+2.*vj(i,j,k)*dt*om*costet
c . . .  Давление через гидростатическую плотность ...
            p1=an6(i+1,j,k)*ro(i+1,j,k)*rg/am_s(an1(i+1,j,k),
     *         an2(i+1,j,k),an3(i+1,j,k))
            p2=an6(i-1,j,k)*ro(i-1,j,k)*rg/am_s(an1(i-1,j,k),
     *         an2(i-1,j,k),an3(i-1,j,k))
c . . .    суммарная плотность для расчета градиента давления
C           ansp=(an1(i+1,j,k)+an2(i+1,j,k)+an3(i+1,j,k))
C           ansl=(an1(i-1,j,k)+an2(i-1,j,k)+an3(i-1,j,k))
C           p1=ansp*bk*an6(i+1,j,k)
C           p2=ansl*bk*an6(i-1,j,k)
C . . .
            f=f-dt*(p1-p2)/(rk*ro(i,j,k)*2.*dtet)
            f=f+vj(i,j,k)**2.*cos_t/(sin_t*rk)*dt
c . . .    ионное трение
            slag=0.89*(an6(i,j,k)+pgi(6,k,i,j))**.37
            caion=(10.8*an2(i,j,k)+12.1*an1(i,j,k)+slag*
     *             an3(i,j,k))*pgi(1,k,i,j)*1.66e-23*1.e-10
            cmion=(12.6*an2(i,j,k)+11.6*an1(i,j,k)+7.8*an3(i,j,k))*
     *             pgl(6,k,i,j)*1.66e-23*1.e-10

            caf=caion*pgi(3,k,i,j)
            f=f+(caf+cmion*vim(k,i,j))/ro(i,j,k)*dt
            b=b+(caion+cmion)/ro(i,j,k)*dt
cc
            pa(k+1)=a/(b-c*pa(k))
            pb(k+1)=(c*pb(k)+f)/(b-c*pa(k))
          end do
          vi1(i,j,n)=pb(n)/(1.-pa(n))
          do l=2,n-1
            k=n-l+1
            vi1(i,j,k)=vi1(i,j,k+1)*pa(k+1)+pb(k+1)
          end do
        end do
      end do
      return
      end
      subroutine vj_r_ham(vj1,vj,vi,vr,vjm,an1,an2,an3,an6,ro,
     *                pgl,pgi,rp,r,eddyco,n,n1,n2,kpars,ins,dt)
      dimension vj(n1,n2,n),vj1(n1,n2,n),vi(n1,n2,n),
     *          vr(n1,n2,n),vjm(n,n1,n2),
     *          an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n),
     *          an6(n1,n2,n),ro(n1,n2,n),
     *          dragVjGSM(n1,n2,n),
     *          pgl(kpars,n,n1,n2),pgi(ins,n,n1,n2),
     *          rp(n),r(n),eddyco(n,n1),
     *          pa(31),pb(31),amu(31)
      data am1,am2,am3,bk,om,pi/53.12e-24,46.51e-24,26.56e-24,
     *     1.38e-16,7.27e-5,3.14159/,re/6.371e08/
      data rg/8.31442e07/
c     . . . средняя масса
        am_s(a1,a2,a3)=(32.*a1+28.*a2+16.*a3)/
     *                 (a1+a2+a3)
c     . . . вычисление шагов по долготе и широте
      dfi=2.*pi/n2
      n1m1=n1-1
      dtet=pi/(n1m1)
      pa(2)=0.
      do i=2,n1m1
c
        tet=dtet*(i-1)
        cos_t=cos(tet)
        sin_t=sin(tet)
        do j=1,n2
          costet=.98*cos_t-.196*sin_t*cos(dfi*(j-1))
          pb(2)=vj(i,j,1)
          do k=1,n
            amu(k)=3.34e-6*(an6(i,j,k)**0.71)
            cp=3.5*(an1(i,j,k)+an2(i,j,k))+2.5*an3(i,j,k)
            cp=cp*bk/(am1*an1(i,j,k)+am2*an2(i,j,k)+
     *         am3*an3(i,j,k))
c       . . . Турбулентная вязкость
            amut=2.*eddyco(k,i)/cp
            amu(k)=amu(k)+amut
          end do
          do k=2,n-1
            rk=r(k)+re
            ot=dt/rp(k)
            del=2.*rp(k)*ro(i,j,k)
            a=dt*(amu(k+1)/rp(k+1)+amu(k)/rp(k))/del
            c=dt*(amu(k)/rp(k)+amu(k-1)/rp(k-1))/del
            vp=0.5*(vr(i,j,k)+abs(vr(i,j,k)))
            vm=0.5*(vr(i,j,k)-abs(vr(i,j,k)))
            vm=vm*ot
            vp=vp*ot
            a=a-vm
            c=c+vp
            b=1.+a+c
c . . .    правая часть
c            f=vj(i,j,k)
c        . . . так надо при данной последовательности расщепления
            f=vj1(i,j,k)
            f=f-2.*vi(i,j,k)*dt*om*costet
            jp=j+1
            jl=j-1
            if(j.eq.1)jl=n2
            if(j.eq.n2)jp=1
c . . .  Давление через гидростатическую плотность ...
            p1=an6(i,jp,k)*ro(i,jp,k)*rg/am_s(an1(i,jp,k),
     *         an2(i,jp,k),an3(i,jp,k))
            p2=an6(i,jl,k)*ro(i,jl,k)*rg/am_s(an1(i,jl,k),
     *         an2(i,jl,k),an3(i,jl,k))
c . . .    суммарная плотность для расчета градиента давления
C            ansp=(an1(i,jp,k)+an2(i,jp,k)+an3(i,jp,k))
C           ansl=(an1(i,jl,k)+an2(i,jl,k)+an3(i,jl,k))
C           p1=ansp*bk*an6(i,jp,k)
C           p2=ansl*bk*an6(i,jl,k)
c
            f=f-dt*(p1-p2)/(ro(i,j,k)*rk*sin_t*2.*dfi)
c . . .   нелинейный член
C            f=f-vi(i,j,k)*vj(i,j,k)*cos_t/(rk*sin_t)*dt
            f=f-vi(i,j,k)*vj1(i,j,k)*cos_t/(rk*sin_t)*dt  ! 18.11.98
c . . .    ионное трение
            slag=0.89*(an6(i,j,k)+pgi(6,k,i,j))**.37
            caion=(10.8*an2(i,j,k)+12.1*an1(i,j,k)+slag*
     *             an3(i,j,k))*pgi(1,k,i,j)*1.66e-23*1.e-10
            cmion=(12.6*an2(i,j,k)+11.6*an1(i,j,k)+7.8*an3(i,j,k))*
     *             pgl(6,k,i,j)*1.66e-23*1.e-10
            caf=caion*pgi(4,k,i,j)

            f=f+(caf+cmion*vjm(k,i,j))/ro(i,j,k)*dt
            b=b+(caion+cmion)/ro(i,j,k)*dt
c
            pa(k+1)=a/(b-c*pa(k))
            pb(k+1)=(c*pb(k)+f)/(b-c*pa(k))
          end do
          vj1(i,j,n)=pb(n)/(1.-pa(n))
          do l=2,n-1
            k=n-l+1
            vj1(i,j,k)=vj1(i,j,k+1)*pa(k+1)+pb(k+1)
c         if(abs(vj1(i,j,k)).gt.1.e4) then
c          print*, ' vj1= ',vj1(i,j,k),vj(i,j,k),i,j,k,pgi(4,k,i,j),
c    *     vjm(k,i,j)
c         end if
          end do
        end do
      end do
      return
      end
c
      subroutine vn_t(vi1,vi,rads,n,n1,n2,dt,key)
c     . . . циклическая прогонка вдоль меридиана
c     . . . для компоненты скорости
c     . . .    nm=n1+n1-2
 
 
      dimension vi(n1,n2,n),vi1(n1,n2,n),rads(n)
      allocatable v_n(:),v(:),a(:),b(:),c(:),f(:)
      
      
      data pi/3.1415926/,re/6.371e8/,par/1./
      nm=2*n1-2

	allocate (v_n(nm),v(nm),a(nm),b(nm),c(nm),f(nm))
	
      ns=n1-1
      n_d=n2/2
      dtet=pi/ns
      ot=dt/dtet
      do k=2,n
        rk=rads(k)+re
        do j=1,n_d
c . . . переход к меридиональному кругу
          !  do i=1,n1
            do i=1,ns  ! 1.11.11
      
			v(i)=vi(i,j,k)
              v_n(i)=vi1(i,j,k)
            end do
          !  do i=n1+1,nm
            do i=n1,nm  ! 1.11.11
              v(i)=-vi(i-n1+1,j+n_d,k)
              if(key.eq.0) then
                v_n(i)=-vi1(i-n1+1,j+n_d,k)
              else
                v_n(i)=vi1(i-n1+1,j+n_d,k)
              end if
            end do
            do i=1,nm
              im=i-1
              ip=i+1
              if(i.eq.nm) ip=1
              if(i.eq.1) im=nm
              vm=(v(i)-abs(v(i)))*.5
              vp=(v(i)+abs(v(i)))*.5
              a(i)=vp*ot*par/rk
              c(i)=-vm*ot*par/rk
              b(i)=1.+a(i)+c(i)
              f(i)=v_n(i)-vm*(1.-par)*ot/rk*(v_n(ip)-v_n(i))-
     *            vp*(1.-par)*ot/rk*(v_n(i)-v_n(im))
            end do
            call cyclp(a,b,c,f,v_n,nm)
            do i=1,n1
              vi1(i,j,k)=v_n(i)
            end do
            do i=n1+1,nm
              if(key.eq.0) then
                vi1(i-n1+1,j+n_d,k)=-v_n(i)
              else
                vi1(i-n1+1,j+n_d,k)=v_n(i)
              end if
            end do
        end do
       end do
       deallocate (v_n,v,a,b,c,f)
       return
       end
c     . . . Расщепление вдоль азимута
c     . . . циклическая прогонка
      subroutine vn_fi(vi1,vj,rads,n,n1,n2,dt)
     
      dimension vi1(n1,n2,n),vj(n1,n2,n),rads(n)
     
      allocatable v_n(:),a(:),b(:),c(:),f(:)
      allocate(v_n(n2),a(n2),b(n2),c(n2),f(n2))
c     . . . par - параметр неявности схемы
       data pi/3.1415926/,re/6.371e8/,par/1./
       ns=n1-1
       dfi=2.*pi/n2
       ot=dt/dfi
       do k=2,n
        rk=rads(k)+re
        do i=2,ns
          tet=pi*(i-1)/ns
          del=rk*sin(tet)
          do j=1,n2
           v_n(j)=vi1(i,j,k)
          end do
          do j=1,n2
            jm=j-1
            jp=j+1
            if(j.eq.1) jm=n2
            if(j.eq.n2) jp=1
            vm=(vj(i,j,k)-abs(vj(i,j,k)))*.5
            vp=(vj(i,j,k)+abs(vj(i,j,k)))*.5
            a(j)=vp*ot*par/del
            c(j)=-vm*ot*par/del
            b(j)=1.+a(j)+c(j)
            f(j)=v_n(j)-vm*(1.-par)*ot/del*(v_n(jp)-v_n(j))-
     *           vp*(1.-par)*ot/del*(v_n(j)-v_n(jm))
          end do
          call cyclp(a,b,c,f,v_n,n2)
          do j=1,n2
            vi1(i,j,k)=v_n(j)
          end do
        end do
       end do
       deallocate(v_n,a,b,c,f)
       return
       end
c
c    . . . производная по времени на 3-х временных слоях
      subroutine vrprim_s(vr,vi,vj,vi1,vj1,ro0,ro,ro1,
     *                   rp,rads,n,n1,n2,dt,nn)
 
      dimension ro1(n1,n2,n),ro(n1,n2,n),ro0(n1,n2,n)
     *          ,rp(n),rads(n),vi1(n1,n2,n),vj1(n1,n2,n)
     *          ,vi(n1,n2,n),vj(n1,n2,n),vr(n1,n2,n)
      allocatable f(:),s(:,:,:),s1(:,:,:),s0(:,:,:),a(:)
      allocate (f(n),s(n1,n2,N),s1(n1,n2,N),s0(n1,n2,N),a(3))
      data pi,re/3.1415926,6.371e 8/,key/1/
      ns=n1-1
      ili=0
      np=n2-1
      nm1=nn-1
      da=pi*2./n2
      df=pi/ns
c     . . . плотность на трех временных слоях
      do 11 i=1,n1
       do 12 j=1,n2
        do 13 k=1,nn
          s(i,j,k)=alog(ro(i,j,k))
          s1(i,j,k)=alog(ro1(i,j,k))
          s0(i,j,k)=alog(ro0(i,j,k))
   13    continue
   12  continue
   11 continue
      do 1 i=2,ns
       tet=df*(i-1)
       sin t=sin(tet)
       sintM1=sin(tet-df)
       sintP1=sin(tet+df)
       cot=1./tan(tet)
       do 2 j=1,n2
          jp=j+1
          jm=j-1
          if(j.eq.n2) jp=1
          if(j.eq.1) jm=n2
!!!!!!!!!!!!!!!!!!!!!!!!!! 30.10.18 !!!!!!!!!!!!!!!
          f(1)=(vr(i,j,2)-vr(i,j,1))/rp(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do 3 k=2,nn
            rk=rads(k)+re
            f(k)=(vr(i,j,k)-vr(i,j,k-1))/rp(k-1)
c    . . . временная производная по s
           if(key.eq.1) then
              f(k)=f(k)+2.*(s1(i,j,k)-s(i,j,k))/dt
              key=0
           else
!              f(k)=f(k)+(s1(i,j,k)-s0(i,j,k))/dt
!            f(k)=f(k)+(s1(i,j,k)-s(i,j,k)+s(i,jm,k)-s0(i,jm,k))/dt !ANALOG KABARE
             f(k)=f(k)+(s1(i,j,k)-s(i,j,k)+s(i,j,k-1)-s0(i,j,k-1))/dt !ANAL OG K.
!             f(k)=f(k)+2.*(s1(i,j,k)-0.5*(s(i,j,k)+s0(i,j,k)))/dt  !12.11.2011
           end if
!           div=.5*(vi1(i+1,j,k)-vi1(i,j,k)+vi(i,j,k)-vi(i-1,j,k))/
!     *          df+(vi1(i,j,k)+vi1(i,j,k))*.5*cot  ! 18.11.98
!c    *          df+vi1(i,j,k)*cot
                !!!!!!!!!!!! new var  27.07.18  !!!!!!!!!
           div=.5*(vi1(i+1,j,k)*sintP1/sint-vi1(i,j,k)+
     *          vi(i,j,k)-vi(i-1,j,k)*sintM1/sint)/df
           !!!!!!!!!!!! new var  !!!!!!!!!         
           div=div+.5/(da*sin t)*(vj1(i,jp,k)-vj1(i,j,k)+
     *         vj(i,j,k)-vj(i,jm,k))
           vipl=(vi1(i,j,k)+abs(vi1(i,j,k)))*.5
           vim=(vi1(i,j,k)-abs(vi1(i,j,k)))*.5
           vjpl=(vj1(i,j,k)+abs(vj1(i,j,k)))*.5
           vjm=(vj1(i,j,k)-abs(vj1(i,j,k)))*.5
           grad=vipl*(s(i,j,k)-s(i-1,j,k))+
     *          vim*(s(i+1,j,k)-s(i,j,k))
           grad=grad/df+(vjpl*(s(i,j,k)-s(i,jm,k))+
     *          vjm*(s(i,jp,k)-s(i,j,k)))/(sin t*da)
           f(k)=(f(k)+(div+grad)*2./rk)*rp(k)
    3   continue
        sig=1.
        if(f(nm1).lt.0.) sig=-1.
        del=(1.+sig)*(s(i,j,nm1)-s(i,j,nn-2))*(rp(nm1)/rp(nn-1))
     *     +(1.-sig)*(s(i,j,nn)-s(i,j,nm1))
        vr(i,j,nm1)=-f(nm1)/del
        vr(i,j,nn)=vr(i,j,nm1)
        do 4 l=2,nm1
         k=nn-l+1
         right=vr(i,j,k+1)+f(k)
         sig=1.
         if(right.lt.0.) sig=-1.
         del=1.-(1.+sig)*(s(i,j,k)-s(i,j,k-1))*(rp(k)/rp(k-1))-
     *       (1.-sig)*(s(i,j,k+1)-s(i,j,k))
!         if(del.lt.1.) then
!             print 119,s(i,j,k+1),s(i,j,k),s(i,j,k-1),i,j,k
!  119        format(' VRPRIM_s  DEL incorrect',1p3e10.2,3i4)
!         end if
         vr(i,j,k)=right/del
         amod=rp(k)/dt
         if(abs(vr(i,j,k+1)).gt.amod) then
          ili=ili+1
          if(ili.ge.1) then
           print 100,vr(i,j,k),s(i,j,k+1),f(k),del,i,j,k
           ili=0
          end if
         end if
    4   continue
!!!!!!!!!!!!!! vr at 1st point !!!!
            k=1
            rk=rads(k)+re
 
c    . . . временная производная по s
              f(k)=f(k)+2.*(s1(i,j,k)-s(i,j,k))/dt
              
           div=.5*(vi1(i+1,j,k)*sintP1/sint-vi1(i,j,k)+
     *          vi(i,j,k)-vi(i-1,j,k)*sintM1/sint)/df
           !!!!!!!!!!!! new var  !!!!!!!!!         
           div=div+.5/(da*sin t)*(vj1(i,jp,k)-vj1(i,j,k)+
     *         vj(i,j,k)-vj(i,jm,k))
           vipl=(vi1(i,j,k)+abs(vi1(i,j,k)))*.5
           vim=(vi1(i,j,k)-abs(vi1(i,j,k)))*.5
           vjpl=(vj1(i,j,k)+abs(vj1(i,j,k)))*.5
           vjm=(vj1(i,j,k)-abs(vj1(i,j,k)))*.5
           grad=vipl*(s(i,j,k)-s(i-1,j,k))+
     *          vim*(s(i+1,j,k)-s(i,j,k))
           grad=grad/df+(vjpl*(s(i,j,k)-s(i,jm,k))+
     *          vjm*(s(i,jp,k)-s(i,j,k)))/(sin t*da)
           f(k)=(f(k)+(div+grad)*2./rk)*rp(k)

         right=vr(i,j,k+1)+f(k)
         
         sig=-1.
         del=1.-(1.-sig)*(s(i,j,k+1)-s(i,j,k))
         vr(i,j,k)=0.5*right/del
!!!!!! 1st point end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!! median smoothing !!!!!!!!!!!!!!!!!!!!!!!!!!
         do k=2,nm1
             a(1:3)=vr(i,j,k-1:k+1)
	     vr(i,j,k)=amed3(a)
	   end do
       
    2  continue
    1 continue
      
  100 format(' VR на данной высоте слишком велико!   vr=',
     *        1pe10.2,' ',3e10.3,3i4)
      deallocate( f,s,s1,s0,a)
       return
       end
c
      function amed3(a)
      dimension a(3)
      amin=a(1)
      amax=a(1)
      imin=1
      imax=1
      do i=2,3
       if(a(i).le.amin) then
          amin=a(i)
          imin=i
       end if

       if(a(i).gt.amax) then

          amax=a(i)
          imax=i
        end if
      end do
      imed=6-(imax+imin)
      if (imed.gt.3) imed=3
      amed3=a(imed)
      return
      end
