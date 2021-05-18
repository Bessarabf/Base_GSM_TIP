!  subroutine sdizkn_bas, gson, kmoldn, o2pro, progon, tn_ic, tn_jc, cyclp
      subroutine sdizkn_bas(an1,an2,an3,an11,an21,an31,an6,vr,
     *                  vi,vj,ro,rp,r,g,n,n1,n2,dt,ctd,ro1,
     *                  solu,gkoor,delta,nsu,dtet,uts,ddolg)
      USE mo_bas_gsm, ONLY:amo2,amn2,amo,qdis
      dimension an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n)
     *         ,ctd(n),an11(n1,n2,n),an21(n1,n2,n),an31(n1,n2,n)
     *         ,vr(n1,n2,n),r(n),rp(n),g(n)
     *         ,ro(n1,n2,n),ro1(n1,n2,n)
     *         ,an6(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n)
     *         ,solu(nsu),gkoor(2,n1,n2)
!      data amo2,amn2,amo/53.12e-24,46.51e-24,26.56e-24/
c
	allocatable q(:,:,:) ! dissosiation sourse
	allocate (q(n1,n2,n))
      n11=n1-1
      lk=17
!!!!!!!!!!!!!!!!!!!!!!! experiment !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loov= 20   ! 20 - base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     . . . Изменена для высокой активности
cc    loov=28  
cc      lov=21 ! зимний вариант
       lov=20  ! основной вариант
c      lov=26
c     lov=n

c    . . .  источник фотодиссоциации (q)
      do i=1,n1
        do j=1,n2
!          call fqsmen(an1,an2,an3,an6,n,n1,i,j,solu,q,
!     *                gkoor,r,delta,n2,nsu,uts)
          do k=1,n
		  q(i,j,k)=qdis(1,i,j,k)
	    end do
        end do
      end do

      n4=n
c         O    con,
      call gson  (an1,an31,an2,an3,vr,an6,q,rp,r,g,lov,dt,
     *            ctd,ro1,n1,n2,n ,vi,vj)
      call boskli(an31,n,n1,n2)
      call barsos(an31,an6,rp,g,amo,n,n1,n2,lov)
c . . . Циклическая прогонка
!! comment 30.07.18 non explicit scheme
       call tn_jc(an31,vj,
     *            r,n,n1,n2,dt,n4)
       call boskli(an31,n,n1,n2)
       call tn_ic(an31,vi,
     *           r,n,n1,n2,dt,n4)
       call bongl(an31,n,n1,n2)
!!      call ficon(an31,r,n,n1,n2,dt,ddolg,vj,dtet,vi)
!!      call boskli(an31,n,n1,n2)
!!  end non explicit scheme 
!!          explicit scheme 
!      call tngojm(an31,vj,
!     *     r,n,n1,n2,dt,n)
!      call boskli(an31,n,n1,n2)
!      call tngoim_a(an31,vi,     ! advection across poles
!     *     r,n,n1,n2,dt,n)
!      call bongl(an31,n,n1,n2)
c . . . Old Var
c      call boskli(an31,n,n1,n2)
c      call barsos(an31,an6,rp,g,amo2,n,n1,n2,lov)
c . . . явная схема переноса через полюс
c       call tetcon(an31,r,n,n1,n2,dt,dtet,vi,vj,n4)
c      call tetcon_a(an31,vi,r,n,n1,n2,dt,n4)
c      call boskli(an31,n,n1,n2)
c
c         O2    con.
c
!       call gsoom (an1,an11,an2,an3,vr,an6,q,rp,r,g,loov,dt,
!     *              ctd,ro1,n1,n2,n ,vi,vj)
c     . . . прогонка
       call o2pro  (an11,an1,an2,an3,an6,vr,vi,vj,
     *               q,ctd,r,rp,g,n1,n2,n,dt)
       call boskli(an11,n,n1,n2)
!  . . .   Корректировка О2 с ',loov,' точки'
       call barsos(an11,an6,rp,g,amo2,n,n1,n2,loov)! loov)
!!   explicit scheme 
!       call tngojm(an11,vj,
!     *     r,n,n1,n2,dt,n)
!       call boskli(an11,n,n1,n2)
!       call tngoim_a(an11,vi,     ! advection across poles
!     *      r,n,n1,n2,dt,n) 
!       call bongl(an11,n,n1,n2)

!       call ficon(an11,r,n,n1,n2,dt,ddolg,vj,dtet,vi)
!       call boskli(an11,n,n1,n2)
c       call tetcon(an11,r,n,n1,n2,dt,dtet,vi,vj,n4)
c      call tetcon_a(an11,vi,r,n,n1,n2,dt,n4)
c      call boskli(an11,n,n1,n2)
c      call barsos(an11,an6,rp,g,amo2,n,n1,n2,loov)
c      call bongl(an11,n,n1,n2)
!! end of explicit scheme
c . . . Циклическая прогонка
       call tn_jc(an11,vj,
     *            r,n,n1,n2,dt,n)  ! loov)
       call boskli(an11,n,n1,n2)
       call tn_ic(an11,vi,
     *            r,n,n1,n2,dt,n)  ! loov)
       call bongl(an11,n,n1,n2)
!       call boskli(an11,n,n1,n2)
c          N2 con.
       key=0
c      lkk=lk
       lkk=n-1
       ot=amo2/amn2
       ot1=amo/amn2
       do 1 i=2,n11
        do 2 j=1,n2
         do 3 k=2,n-1
          tot=(ro1(i,j,k))/amn2
          sum=ot*an11(i,j,k)+ot1*an31(i,j,k)
          an21(i,j,k)=tot-sum
c         . . . Автоматическая проверка
          prov=0.5*an31(i,j,k)
c          prov=0.3*an31(i,j,k)
          if(an21(i,j,k).lt.prov) then
           if(k.lt.lkk) then
            key=1
            lkk= k-1
            if(lkk.LE.2) lkk=2 
            ii=i
            jj=j
           end if
          end if
    3    continue
    2   continue
    1  continue
       if (key.eq.1)
     *     print *,' N2 subtract until',lkk, ' point',ii,jj
c    *     print 100,an21(ii,jj,lkk),ii,jj,lkk
      call barsos   (an21,an6,rp,g,amn2,n,n1,n2,lkk)
c     call nts(an21,n,n1,n2,n4)
      call boskli(an21,n,n1,n2)
  100 format(' ОТРИЦАТЕЛЬНАЯ [N2]=',e10.2,3i4)
      
	deallocate (q)

      return
      end
!------------------------------------------------------------------------------
      subroutine gson(an1,an31,an2,an3,vr,an6,q,
     *                  rp,r,g,lv,dt,ctd,ro,n1,n2,n,vi,vj)
c . . .
      dimension an1(n1,n2,n),an2(n1,n2,n),
     *          an3(n1,n2,n),an6(n1,n2,n),an31(n1,n2,n),
     *          vr(n1,n2,n),q(n1,n2,n),rp(n),r(n),g(n),ctd(n),
     *          ro(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n),
     *          dp(3,3)
      allocatable aal(:),bbc(:),dd31(:),dd32(:),ql(:),
     *          fp(:),ap(:),bp(:),ep(:),ed(:),b2(:),
     *          b3(:),cp(:),dpp(:),pfu(:),ppo(:)
      data am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
      data bk/1.38e-16/,pi/3.14159/,re/6.371e08/
      
	
      allocate (aal(n),bbc(n),dd31(n),dd32(n),ql(n),
     *          fp(n),ap(n),bp(n),ep(n),ed(n),b2(n),
     *          b3(n),cp(n),dpp(n),pfu(n),ppo(n))


      if(n.lt.lv) print *,'mistake on gsozkg. n can''t low lv'
      np=lv-1
      n11=n1-1
      n21=n2
      dtet=pi/(n1-1)
      dfi=2.*pi/n2
      gam=1.e-20
c*******
      do 40 i=2,n11
       tetl=dtet*(i-1)
       do 401 j=1,n21
        jp=j+1
        jl=j-1
        if(j.eq.1)jl=n2
        if(j.eq.n21)jp=1
c
        do 1 k=1,lv
            ams=(32.*an1(i,j,k)+an2(i,j,k)*28.+16.*an3(i,j,k))
     *         /(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))*1.67e-24
            cono2=an1(i,j,k)
            conn2=an2(i,j,k)
            cono=an3(i,j,k)
            tempp=an6(i,j,k)
            rog=ro(i,j,k)
c
c           call kmold(dp,cono2,conn2,cono,tempp)
            call kmoldn(dp,cono2,conn2,cono,tempp,rog)
            hh=bk*tempp/g(k)
            hcp=hh/ams
            ho=hh/am3
            hoo=hh/am1
            hnn=hh/am2
            kv=k+1
            kn=k-1
            if(k.eq.lv)kv=lv
            if(k.eq.1)kn=1
            gtem=(an6(i,j,kv)-an6(i,j,kn))/(r(kv)-r(kn))
            aal(k)=ctd(k)+dp(3,3)
            bbc(k)=ctd(k)*(1./hcp+1./tempp*gtem)+
     *      dp(3,3)*(1./ho+1./tempp*gtem)
c           bbc(k)=bbc(k)-vr(i,j,k)*0.
cc          b2(k)=dp(3,1)*(1./hoo+1./tempp*gtem)
cc            b3(k)=dp(3,2)*(1./hnn+1./tempp*gtem)
cc            b2(k)=b2(k)/an3(i,j,k)*an1(i,j,k)
cc            b3(k)=b3(k)/an3(i,j,k)*an2(i,j,k)
cc            goo=(an1(i,j,kv)-an1(i,j,kn))/(r(kv)-r(kn))
cc   *        *    dp(3,1)/an3(i,j,k)
cc            gnn=(an2(i,j,kv)-an2(i,j,kn))/(r(kv)-r(kn))
cc   *        *    dp(3,2)/an3(i,j,k)
cc            b2(k)=b2(k)+goo
cc            b3(k)=b3(k)+gnn
              rt=an6(i,j,k)*(r(kv)-r(kn))
              ht=(an1(i,j,kv)*an6(i,j,kv)-
     *            an1(i,j,kn)*an6(i,j,kn))/rt
              hb=an1(i,j,k)/hoo
                hs=ht+hb
c               hm=hs/hb
                hm=hs/ht
cc            if(abs(hm).le.0.1)hs=0.
                b2(k)=dp(3,1)*hs/an3(i,j,k)

              ht=(an2(i,j,kv)*an6(i,j,kv)-
     *            an2(i,j,kn)*an6(i,j,kn))/rt
              hb=an2(i,j,k)/hnn
                hs=ht+hb
c               hm=hs/hb
                hm=hs/ht
cc            if(abs(hm).le.0.1)hs=0.
                b3(k)=dp(3,2)*hs/an3(i,j,k)
              dd31(k)=dp(3,1)
              dd32(k)=dp(3,2)
 1        continue
c*******
          do 3 k=2,np
            alf=9.9e-34*exp(470./an6(i,j,k))
            bet=1.1e-34*exp(510./an6(i,j,k))
            pp=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
            if(q(i,j,k).ne.0.)bet=0.
            qq=2.*q(i,j,k)*an1(i,j,k)
            ql(k)=2.*alf*an3(i,j,k)*pp+bet*an1(i,j,k)*pp
     *      +2.*gam*an3(i,j,k)
c*******
            aa=0.5/dfi/(r(k)+re)/sin(tetl)
            dd=aa*(vj(i,jp,k)-vj(i,jl,k))
            aa=0.5/dtet/(r(k)+re)/sin(tetl)
c            dd=dd+aa*(vi(i+1,j,k)*sin(tetl+dtet)
c     *         -vi(i-1,j,k)*sin(tetl-dtet))
c . . .  Дивергенция V с котангенсом:
            dd=dd+vi(i,j,k)/(r(k)+re)*(cos(tetl)/sin(tetl))+
     *         (vi(i+1,j,k)-vi(i-1,j,k))/(r(k)+re)/dtet*0.5
            alf=dd
            fp(k)=qq
            if(alf.gt.0.)ql(k)=ql(k)+alf
            if(alf.le.0.)fp(k)=fp(k)-alf*an3(i,j,k)
 3        continue
c*******
          fp(1)=0.
          fp(lv)=0.
          do 2 k=1,np
            dx=r(k+1)-r(k)
            bbc2=(bbc(k)+bbc(k+1))/2.
            aal2=(aal(k)+aal(k+1))/2.
            aaa=(bbc2/aal2+bbc(k+1)/aal(k+1))/4.*dx
            bbb=(bbc2/aal2+bbc(k)/aal(k))/4.*dx
            if (aaa.gt.40.) then
             print 777, i,j,k,aaa
             aaa=40.
            endif
            ed(k)=exp(aaa)*aal2/dx
            if (bbb.gt.40.) then
             print 778, i,j,k,bbb
             bbb=40.
            endif
            ep(k)=exp(-bbb)*aal2/dx
            ppm=(vr(i,j,k+1)-abs(vr(i,j,k+1)))/2.
            ppp=(vr(i,j,k)+abs(vr(i,j,k)))/2.
            if((vr(i,j,k+1).lt.0.).and.(vr(i,j,k).gt.0.))ppm=0.
            ppma=(b2(k+1)-abs(b2(k+1)))/2.
            pppa=(b2(k)+abs(b2(k)))/2.
            ppmb=(b3(k+1)-abs(b3(k+1)))/2.
            pppb=(b3(k)+abs(b3(k)))/2.
            ed(k)=ed(k)-ppm-ppma-ppmb
            ep(k)=ep(k)+ppp+pppa+pppb
 2        continue
c*******
          ap(1)=0.
          bp(1)=an3(i,j,1)
          zz=ed(1)
c*******
          do 6 k=2,np
            ddx=(r(k+1)-r(k-1))/2.
            dtt=dt
            cp(k)=ddx*(1./dtt+ql(k))
            fp(k)=(an3(i,j,k)/dtt+fp(k))*ddx
            aaa=ep(k)+cp(k)+zz
            zzz=ed(k)*(cp(k)+zz)/aaa
            ap(k)=ed(k)/aaa
            bp(k)=(fp(k)+ep(k-1)*bp(k-1))/aaa
            zz=zzz
 6        continue
c*******
          hgr=r(lv)-r(np)
          ho=bk*an6(i,j,lv)/g(lv)/am3
          ppo(np)=bp(np)/(1.-ep(np)/ed(np)*ap(np))
          hod=1./ho
          ppo(lv)=ppo(np)*exp(-hgr*hod)
          an31(i,j,lv)=ppo(lv)
c*******
          do 7 jj=2,lv
            k=lv-jj+1
            ppo(k)=ap(k)*ppo(k+1)+bp(k)
            an31(i,j,k)=(ppo(k))
 7        continue
401   continue
 40   continue
  777 format(' GSO. MISTAKE  aaa ',3i4,1pe10.2)
  778 format(' GSO. MISTAKE  bbb ',3i4,1pe10.2)

      deallocate (aal,bbc,dd31,dd32,ql,
     *          fp,ap,bp,ep,ed,b2,
     *          b3,cp,dpp,pfu,ppo)


      return
      end
!------------------------------------------------------------------------------
      subroutine kmoldn(dp,an1,an2,an3,an6,ro)
c
c     . . p/p ras. k. mol. dif. to trexkomponent.
c     sredi   v odnoj tocke
      dimension
     *          dp(3,3)
      data am1,am2,am3 /53.12e-24,46.51e-24,26.56e-24/,
     *     bk,pi/1.38e-16,3.1415296/,
     *     d1,d2,d3/3.6e-8,3.7e-8,2.e-8/
      coef=3./8.*sqrt(0.5*bk/pi)
      sum=an1+an2+an3
c     . . . otnosheie mol mass
      otm12=am1/am2
      otm13=am1/am3
      otm23=am2/am3
c     . . . srednij diametr vzaim.
      sig12=0.5*(d1+d2)
      sig13=0.5*(d1+d3)
      sig23=0.5*(d2+d3)
      sig12=sig12*sig12
      sig13=sig13*sig13
      sig23=sig23*sig23
c     . . . obrat. massa
      obr1=1./am1
      obr2=1./am2
      obr3=1./am3
      b12=sqrt(obr1+obr2)
      b13=sqrt(obr1+obr3)
      b23=sqrt(obr2+obr3)
      sqt=sqrt(an6)
      value=coef*sqt
      d12=value/sig12*(b12/sum)
      d13=value/sig13*(b13/sum)
      d23=value/sig23*(b23/sum)
      d32=d23
      ro1=an1*am1
      ro2=an2*am2
      ro3=an3*am3
      znam=(an1*d23+an2*d13+an3*d12)*ro
      rez=sum/znam
      dp(1,1)=d13*d12*(ro2+ro3)*rez
      dp(1,2)=-d23*d12*ro1*rez/otm12
      dp(1,3)=-d23*d13*ro1*rez/otm13
      dp(2,1)=-d12*d13*ro2*rez*otm12
      dp(2,2)=d23*d12*(ro1+ro3)*rez
      dp(2,3)=-d23*d13*ro2*rez/otm23
      dp(3,1)=-d12*d13*ro3*rez*otm13
      dp(3,2)=-d23*d12*ro3*rez*otm23
      dp(3,3)=d23*d13*(ro2+ro1)*rez
      return
      end
!------------------------------------------------------------------------------
! new ver 01_11_18
      subroutine o2pro (an11,an1,an2,an3,an6,vr,vi,vj,
     *                  q,ctd,rads,rp,g,n1,n2,n,dt)
c     O2 в приближении малой компоненты  (прогонка)

      dimension an1(n1,n2,n),an2(n1,n2,n),
     *          an3(n1,n2,n),an6(n1,n2,n),an11(n1,n2,n),
     *          vr(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n),
     *          q(n1,n2,n),rads(n),rp(n),g(n),ctd(n)
     *          
      allocatable a(:),b(:),c(:),f(:),cmd(:)
     *         ,h(:),alf(:),bet(:),hsr(:),c o2(:)
      data am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
     *     bk/1.38e-16/,gam/1.e-20/,pi/3.141592/,re/6.371e8/
	allocate (a(n),b(n),c(n),f(n),cmd(n)
     *         ,h(n),alf(n),bet(n),hsr(n),c o2(n))
      const=bk/am1
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
         c o2(1)=an1(i,j,1)
         do k=1,n
          h(k)=const*an6(i,j,k)/g(k)
          sum=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
          ams=(am1*an1(i,j,k)+am2*an2(i,j,k)+am3*an3(i,j,k))/sum
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
 !         d12=0.829e17/sum1*an6(i,j,k)**0.724*epok
 !         d13=0.969e17/sum2*an6(i,j,k)**0.774*epok
 !         obr=(an2(i,j,k)/d12+an3(i,j,k)/d13)/(an2(i,j,k)+an3(i,j,k))
 !         cmd(k)=1./obr
          cmd(k)=3.e17/sum*sqrt(an6(i,j,k))  ! к-т м.диффузии NO
          alf(k)=cmd(k)/(cmd(k)+ctd(k))
          bet(k)=ctd(k)/(cmd(k)+ctd(k))
         end do
         c o2(1)=an1(i,j,1)
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
          cmdp=(cmd(k+1)+cmd(k)+ctd(k+1)+ctd(k))*.5
          cmdm=(cmd(k)+cmd(k-1)+ctd(k)+ctd(k-1))*.5
c          vrp=(vr(i,j,k)+abs(vr(i,j,k)))*.5
c          vrm=(vr(i,j,k)-abs(vr(i,j,k)))*.5
c          a(k)=(cmdp/pro-vrm)/rp(k)
c          c(k)=(cmdm/pro+vrp)/rp(k-1)
          a(k)=(cmdp/pro)/rp(k)
          c(k)=(cmdm/pro)/rp(k-1)
          b(k)=a(k)+c(k)+1./dt+p
          clp=0.5*(alf(k+1)/h(k+1)+dtnp/rp(k)+bet(k+1)/hsr(k+1))
          clm=-0.5*(alf(k-1)/h(k-1)+dtnm/rp(k-1)+bet(k-1)/hsr(k-1))
          aprim=(cmd(k+1)+ctd(k+1))/pro*clp
          cprim=(cmd(k-1)+ctd(k-1))/pro*clm
          cl=0.5*(-dtnp/rp(k)+dtnm/rp(k-1))
!          a(k)=a(k)+aprim-vr(i,j,k+1)/pro*.5
!          c(k)=c(k)+cprim+vr(i,j,k-1)/pro*.5 ! 31_10_18
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
          b(k)=b(k)+cl*(cmd(k)+ctd(k))/pro
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
        f(n)=exp(-rp(n-1)/h(n))
        kiss=2
        call progon (c o2,a,b,c,f,n,kiss)
c     . . .
! . . . reASSIGN 

        an11(i,j,1)=c o2(1)
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! ver 20.09.19  left and right  first derivative in 1 point is equal
       subroutine tn_ic(an6,vi,rads,n,n1,n2,dt,n0)
c     . . . циклическая прогонка вдоль меридиана
c     . . .    nm=n1+n1-2


       dimension an6(n1,n2,n),vi(n1,n2,n),rads(n)
       allocatable tn(:),v(:),a(:)
     *,            b(:),c(:),f(:)
       data pi/3.1415926/,re/6.371e8/,par/1./
       nm=n1+n1-2
       allocate (tn(nm),v(nm),a(nm)
     *,          b(nm),c(nm),f(nm))

       ns=n1-1
       n_d=n2/2
       dtet=pi/ns
       ot=dt/dtet
       do k=2,n0
        rk=rads(k)+re
        do j=1,n_d
c . . . переход к меридиональному кругу
            do i=1,n1
              v(i)=vi(i,j,k)
              tn(i)=an6(i,j,k)
            end do
            do i=n1+1,nm
              v(i)=-vi(i-n1+1,j+n_d,k)
              tn(i)=an6(i-n1+1,j+n_d,k)
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
              f(i)=tn(i)-vm*(1.-par)*ot/rk*(tn(ip)-tn(i))-
     *                   vp*(1.-par)*ot/rk*(tn(i)-tn(im))
            end do
            call cyclp(a,b,c,f,tn,nm)
            do i=1,n1
              an6(i,j,k)=tn(i)
            end do
            do i=n1+1,nm
              an6(i-n1+1,j+n_d,k)=tn(i)
            end do
        end do
       end do
       deallocate (tn,v,a
     *,          b,c,f)
       return
       end
!------------------------------------------------------------------------------
       subroutine tn_jc(an6,vj,rads,n,n1,n2,dt,n0)
c     . . . циклическая прогонка
  
       dimension an6(n1,n2,n),vj(n1,n2,n),rads(n)
       allocatable tn(:),a(:),b(:),c(:),f(:)
       allocate (tn(n2),a(n2),b(n2),c(n2),f(n2)) 
       data pi/3.1415926/,re/6.371e8/,par/1./
       nm=n2
       ns=n1-1
       dfi=2.*pi/n2
       ot=dt/dfi
       do k=2,n0
        rk=rads(k)+re
        do i=2,ns
          tet=pi*(i-1)/ns
          del=rk*sin(tet)
          do j=1,n2
           tn(j)=an6(i,j,k)
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
            f(j)=tn(j)-vm*(1.-par)*ot/del*(tn(jp)-tn(j))-
     *                 vp*(1.-par)*ot/del*(tn(j)-tn(jm))
          end do
          call cyclp(a,b,c,f,tn,n2)
          do j=2,n2
            an6(i,j,k)=tn(j)
          end do
!!!     left and right  first derivative in 1 point is equal
          an6(i,1,k)=(tn(2)+tn(n2))*.5
        end do
       end do
       deallocate (tn,a,b,c,f)
       return
       end
!------------------------------------------------------------------------------
c . . . циклическая прогонка
c . . . a, b, c, f - к-ты трехточки
c . . . a(i)y(i-1)-b(i)y(i)+c(i)y(i+1)=-f(i)
      subroutine cyclp(a,b,c,f,y,nx)
c . . . nx1=nx+1
      dimension a(nx),b(nx),c(nx),f(nx),y(nx)
      allocatable alf(:),bet(:),gam(:)
     *         ,p(:),q(:)
      allocate (alf(nx+1),bet(nx+1),gam(nx+1)
     *         ,p(nx+1),q(nx+1))
       alf(2)=c(2)/b(1)
       gam(2)=a(1)/b(1)
       bet(2)=f(1)/b(1)
       do i=2,nx
         ab=b(i)-a(i)*alf(i)
         alf(i+1)=c(i)/ab
         gam(i+1)=a(i)*gam(i)/ab
         bet(i+1)=(a(i)*bet(i)+f(i))/ab
       end do
       p(nx-1)=bet(nx)
       q(nx-1)=alf(nx)+gam(nx)
       do j=nx-2,1,-1
         p(j)=alf(j+1)*p(j+1)+bet(j+1)
         q(j)=alf(j+1)*q(j+1)+gam(j+1)
       end do
       del=1-alf(nx+1)*q(1)-gam(nx+1)
       y(nx)=(bet(nx+1)+alf(nx+1)*p(1))/del
       do i=1,nx-1
        y(i)=p(i)+y(nx)*q(i)
       end do
       deallocate(alf,bet,gam,p,q)
       return
       end
!------------------------------------------------------------------------------