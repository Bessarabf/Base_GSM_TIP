C Temperature calculation
C 10.05.2018 - Joul source send to HAMMONIA
C list of programm:
C - HEATPO
C - TNPOT
C - IJOULP
C - ANOIK
C - COLD
C - DJOULP
C - FREAC
C - E

c    . . . cicle_prog alog i и j
      subroutine heatpo_bas(pgl,pgi,parj,solet,solu,nsu,nse,
     *           kpars,rads,g,nh,gkoor,its,ddolg,dtet,
     *           mass,delta,day,uts,
     *           tau,dts,ctd,vim,vid,vir,ids,ins,
     *           an1,an2,an3,anco2,an6,an61,ro,vi,vj,vr)
      dimension pgl(kpars,nh,its,ids)
     *,         pgi(ins,nh,its,ids),parj(nh,its,ids)
     *,         solet(nse),solu(nsu),rads(nh),mass(30),ctd(nh)
     *,         vim(nh,its,ids),vid(nh,its,ids)
     *,         gkoor(2,its,ids),vir(nh,its,ids)
     *,         an1(its,ids,nh),an2(its,ids,nh)
     *,         an3(its,ids,nh),an6(its,ids,nh),ro(its,ids,nh)
     *,         an61(its,ids,nh),vi(its,ids,nh),vj(its,ids,nh)
     *,         vr(its,ids,nh),g(nh),anco2(its,ids,nh)
     *,         apm(7)
      integer day
      data apm/7*11./,fa,fs/130.,130./
      data pi,om/3.14159,7.27e-5/
!      if(ids.gt.24) then
!       print *,'HEATPO: ids is too small '
!       stop
!      end if
cc    . . .  TN on 526 км according MSIS-90
c     if(mass(10).eq.30) then
c       do   i =  1,its
c         do   j = 1,ids
c           fig=gkoor(1,i,j)
c           fig=90.-fig
c           dgeo=gkoor(2,i,j)/180.*pi
c           dol=gkoor(2,i,j)
c           tau1=uts+dgeo/om
c           tau1=tau1/3600.
c           if(tau1.gt.24.) tau1=tau1-24.
c           td=day
c           iyd=80*1000+td
c            vis=rads(nh)/1.e5
c           vis=526.
c           call gtd6(iyd,uts,vis,fig,dol,tau1,fa,fs,
c    *                apm,48,dm,tm)
c           pgl(7,nh,i,j)=tm(2)
c           end do
c          end do
c     else
c         print*,' Значения Tn на в.границе только при mass(10)=30'
c         stop
cc    end if
!
!    
      call tnpot_bas(pgl,pgi,an1,an2,an3,an6,an61,vi,vj,vr,
     *              anco2,ro,vim,vid,vir,rads,g,gkoor,ctd,solu,nsu,
     *              kpars,ins,nh,its,ids,delta,uts,dts,mass)
      
c  . . . циклическая прогонка вдоль  fi & tet
C     call tn_jc(an61,vj,
C    *     rads,nh,its,ids,dts,mass(10))
C     call boskli(an61,nh,its,ids)
C     call tn_ic(an61,vi,
C    *     rads,nh,its,ids,dts,mass(10))
!!      call bongl(an61,nh,its,ids)
C     call boskli(an61,nh,its,ids)
C     call tnalt(an61,nh,its,ids,mass(10))
C
c  . . . долготная прогонка с учетом теплопроводности
       call tnl_jc(an61,an1,an2,an3,vj,ctd,
     *             rads,nh,its,ids,dts,mass(10))
       call boskli(an61,nh,its,ids)
c  . . . широтная прогонка
       call tnl_ic(an61,an1,an2,an3,vi,ctd,
     *             rads,nh,its,ids,dts,mass(10))
       call boskli(an61,nh,its,ids)
       call tnalt(an61,nh,its,ids,mass(10))
C
c    . . . Старый вариант
c     call tngojm(an61,vj,
c    *     rads,nh,its,ids,dts,mass(10))
c     call boskli(an61,nh,its,ids)
c     call tngoim_a(an61,vi,             ! Явная схема переноса
c    *     rads,nh,its,ids,dts,mass(10)) ! через полюс
c     call boskli(an61,nh,its,ids)
c     call tnalt(an61,nh,its,ids,mass(10))
      return
      end
c
      subroutine tnpot_bas(pgl,pgi,an1,an2,an3,an6,an61,vi,vj,vr,
     *                 anco2,ro,vim,vid,vir,rads,g,gkoor,ctd,solu,nsu,
     *                 kpars,ins,nh,its,ids,dl,uts,dts,mass)
!      USE mo_ham_gsm, ONLY:qJGSM
      dimension pgl(kpars,nh,its,ids),an1(its,ids,nh),
     *          an2(its,ids,nh),an3(its,ids,nh),anco2(its,ids,nh),
     *          an6(its,ids,nh),an61(its,ids,nh),
     *          vj(its,ids,nh),vi(its,ids,nh),
     *          vr(its,ids,nh),rads(nh),mass(30)
     *,         pgi(ins,nh,its,ids)
     *,         solu(nsu),ctd(nh),g(nh)
     *,         vim(nh,its,ids),vid(nh,its,ids),gkoor(2,its,ids)
     *,         vir(nh,its,ids),alyam(3),ro(its,ids,nh)
     
      allocatable pa(:),pb(:),q(:),qdj(:)
      allocate (pa(NH+5),pb(NH+5),q(NH),qdj(NH))
      data re/6.371e8/,pi/3.1415926/,bk/1.38e-16/,
     *    am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
     *,   om/7.272205e-5/
      cr=pi/180.
      klik=0
      n0=mass(10)
      n0m=n0-1
      itsm1=its-1
      itsm2=its-2
      dteta=pi/itsm1
      ddol=2.*pi/ids
c
c
      do 1 i=2,itsm1
       teta=dteta*(i-1)
       sin t=sin(teta)
       sin p=sin(teta+dteta)
       sin m=sin(teta-dteta)
       do 2 j=1,ids
        an60=pgl(7,1,i,j)
        call ijoulp(q,qdj,pgl,pgi,ctd,rads,g,an60,solu,
     *            gkoor,kpars,ins,nh,its,ids,nsu,i,j,uts,dl,
     *            an1,an2,an3,an6,anco2,vr,vi,vj,vim,vid,vir)
        do 29 k=1,n0m
         rc=(2.5*(an1(i,j,k)+an2(i,j,k))+1.5*an3(i,j,k))*bk
 !        qJGSM(i,j,k)=qdj(k)/rc
         q(k)=q(k)+qdj(k)
   29   continue
c
      pa(2)=0.
      pb(2)=an60
      jm=j-1
      jp=j+1
      if(j.eq.1) jm=ids
      if(j.eq.ids) jp=1
      do 3 k=2,n0m
c
        do 20 k1=1,3
         kv=k-2+k1
          al=18.6*an6(i,j,kv)**0.84*an1(i,j,kv)
          al=al+27.2*an6(i,j,kv)**0.8*an2(i,j,kv)
          al=al+67.1*an6(i,j,kv)**0.71*an3(i,j,kv)
          alyam(k1)=al/(an1(i,j,kv)+an2(i,j,kv)+
     *              an3(i,j,kv))
c         . . . Учет турбулентной теплопроводности
          rocp=(3.5*(an1(i,j,kv)+an2(i,j,kv))+2.5*
     *               an3(i,j,kv))*bk
          alyam(k1)=alyam(k1)+rocp*ctd(kv)
c
   20   continue
c
        anu=3.34e-6*an6(i,j,k)**0.71
c        sum=(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))*bk
        sum1=(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))
        ams=(an1(i,j,k)*am1+an2(i,j,k)*am2 +
     *       an3(i,j,k)*am3)/sum1
        rc=(2.5*(an1(i,j,k)+an2(i,j,k))+1.5*an3(i,j,k))*bk
        abs vr= abs(vr(i,j,k))
        con1=vr(i,j,k)-abs vr
        con2=vr(i,j,k)+abs vr
c
        drad=rads(k+1)-rads(k)
        drad1=rads(k)-rads(k-1)
c
        rk=rads(k)+re
c
        a=((alyam(3)+alyam(2))/drad-rc*con1)/(rads(k+1)-
     *     rads(k-1))
        c=((alyam(2)+alyam(1))/drad1+rc*con2)/(rads(k+1)-
     *     rads(k-1))
        b=a+c+rc/dts
c
        davl0=(vr(i,j,k+1)-vr(i,j,k-1))/(drad+drad1)
c        davl1=(vi(i+1,j,k)*sin p-
c     *         vi(i-1,j,k)*sin m)/dteta*.5
        davl1=(vi(i+1,j,k)-
     *         vi(i-1,j,k))/dteta*.5+vi(i,j,k)*cos(teta)/sin t
        davl2=(vj(i,jp,k)-vj(i,jm,k))/ddol*.5
c        davl=davl0+(davl1+davl2)/(rk*sin t)
        davl=davl0+davl1/rk+davl2/(rk*sin t)
        ha=rc*pgl(7,k,i,j)/dts
c . . .  Вязкий нагрев
        vqi=(vi(i,j,k+1)-vi(i,j,k))/drad
        vqj=(vj(i,j,k+1)-vj(i,j,k))/drad
        vq=(vqi*vqi+vqj*vqj)*anu
c       fk=q(k)+  ha -davl*sum*an6(i,j,k)+vq
         fk=q(k) + ha+vq
c . . .  Учет знака давления
        if(davl.lt.0.) then
          fk=fk-davl*ro(i,j,k)*bk/ams*an6(i,j,k)
        else
          b=b+davl*ro(i,j,k)*bk/ams
        end if
c . . .  Прямая прогонка
        pa(k+1)=a/(b-pa(k)*c)
        pb(k+1)=(c*pb(k)+fk)/(b-pa(k)*c)
    3  continue
c*
       if(mass(22).ne.0) then
c    . . . Учет нелокального нагрева
         q_ot=pgl(16,n0,i,j)+pgl(16,n0,its-i+1,j)
         if(pgl(16,n0,i,j).ge.pgl(16,n0,its-i+1,j)) then
            q_max=pgl(16,n0,i,j)
         else
            q_max=pgl(16,n0,its-i+1,j)
         end if
         q_ot=q_ot/q_max
       ! pot_m=0.05
	   pot_m=0.02                    !   Поток эрг/см2*c-1
c  night flux on upper boundary erg/cm2*c-1 
         if (klik.eq.0) then
           print*,'temp -nonlocal heating. Flux=',pot_m,' erg/cm2*c-1'
           klik=1
         end if
         pot=pot_m*.5*q_ot
         pot=pot/alyam(3)
         rp=rads(n0)-rads(n0-1)
         an61(i,j,n0)=(pb(n0)+pot*rp)/(1.-pa(n0))
       else
         an61(i,j,n0)=pb(n0)/(1.-pa(n0))
       end if
       do 4 l=2,n0m
        k=n0-l+1
        an61(i,j,k)=pa(k+1)*an61(i,j,k+1)+pb(k+1)
        if(an61(i,j,k).lt.100.or.an61(i,j,k).gt.1600.) then
            print 100,an61(i,j,k),k,i,j
        end if
    4  continue
    2 continue
    1 continue
  100 format('GSMTIP: warning Tn=',f6.0,' k=',i3,' i=',i2,' j=',i2)
      deallocate (pa,pb,q,qdj)
      return
      end
!
      subroutine ijoulp(q,qdj,pgl,pgi,ctd,rads,g,an60,solu,
     *                  gkoor,kpars,ins,nh,its,ids,nsu,i,j,uts,
     *                  del,an1,an2,an3,an6,anco2,vr,vi,vj,vim,vid,vir)
      dimension q(nh),g(nh),pgl(kpars,nh,its,ids),gkoor(2,its,ids),
     *          rads(nh),solu(nsu),ctd(nh),sni(6)
     *         ,an1(its,ids,nh),an2(its,ids,nh),an3(its,ids,nh)
     *         ,an6(its,ids,nh),anco2(its,ids,nh),vr(its,ids,nh)
     *         ,vi(its,ids,nh),vj(its,ids,nh)
     *         ,pgi(ins,nh,its,ids)
     *         ,vim(nh,its,ids),vid(nh,its,ids)
     *         ,vir(nh,its,ids),qdj(nh)
      data pi/3.14159265359d0/,om/7.272205e-5/,
c   . . . эффективность для зимы
c    *     r0,r00/3.48,2.94/
c   . . . эффективность=0.6   равноденствие
c    *     r1,r2/3.48,2.94/
c   ..... efficency = 0.3
c     *     r1,r2/1.74,1.47/
c   ..... efficency = 0.45        ! Уменьшена эффективность
!     *     r1,r2/2.61,2.205/      ! для лета
c   ..... efficency = 0.5         ! Уменьшена эффективность
c    *     r1,r2/3.21,2.5 /       ! для лета
c   . . . эффективность=0.4
     *    r1,r2/2.32,1.96/
c   . . . эффективность=0.35
c    *     r1,r2/2.00,1.71/
     *    ,r0,r00/3.48,2.94/      ! для зимы
c   . . . UV - efficency          ! Эффективность для UV
c    *    ,e_dis0,e_dis00/0.3,1.0/! основной вариант
     *    ,e_dis0,e_dis00/0.3,0.8/
      n0=nh-1
c*
      cr=pi/180.
      g shir=gkoor(1,i,j)*cr
      g dol=gkoor(2,i,j)*cr
      g shir=pi/2.-g shir
c*
      f=sin(g shir)*sin(del)+cos(gshir)*cos(del)*
     *  cos(om*(uts-43200.)+g dol)
      hi=acos(f)
      key=1
      q(1)=0.
      tem0=0.
      do 1 k=2,n0
        ano=an1(i,j,k)
        and=an2(i,j,k)
        antr=an3(i,j,k)
        con no=pgl(4,k,i,j)
        tem=an6(i,j,k)
        viam=pgi(3,k,i,j)
        viad=pgi(4,k,i,j)
        viar=pgi(2,k,i,j)
        vimm=vim(k,i,j)
        vimd=vid(k,i,j)
        vimr=vir(k,i,j)
        vin=vi(i,j,k)
        vjn=vj(i,j,k)
        vrn=vr(i,j,k)
        cona=pgi(1,k,i,j)
        conm=pgl(6,k,i,j)
        ti=pgl(8,k,i,j)
        te=pgl(9,k,i,j)
        conco2=anco2(i,j,k)
        
        di=dis mod(ano,tem,g(k),rads(k),solu,nsu,hi,key)*ano
c
        e_dis=0.3          ! Равномерный 
        if(k.le.23) then   ! нагрев
           r=r1          !
        else             !
           r=r2          !
       end if            !
c . . .  для корректировки эффективности:
c . . .  зимой возрастает до 0.6 (только для ЗИМНЕГО солнцестояния!!!)
c . . . Увеличение эффективности в зимнем полушарии
!!        g sh=gkoor(1,i,j)
!!        if(g sh.gt.90.) then  ! Летом без изменения
!!         e_dis=e_dis0
!!         if(k.le.23) then
!!            r=r1
!!         else
!!            r=r2
!!         end if
!!        else
!!          e_dis=e_dis00+(e_dis0-e_dis00)*(g sh)/90.  ! Зима
!!         if(k.le.23) then
!!            r=r0+(r1-r0)*(g sh)/90.  ! Зима
!c           r=r1+(r0-r1)*(g sh-90.)/90.   ! Лето
!!         else
!!            r=r00+(r2-r00)*(g sh)/90.  ! Зима
!c           r=r2+(r00-r2)*(g sh-90.)/90.  ! Лето
!!         end if
!!        end if
        qi=r*(pgl(13,k,i,j)+pgl(14,k,i,j)+
     *        pgl(15,k,i,j)+pgl(16,k,i,j))*1.e-11
        co=cold(ano,and,antr,tem,conco2,rads(k),g(k),ctd(k))
        cik53=anoik(con no,antr,tem)
        che=chem(ano,and,antr,tem)
c*
        pol=(tem+pgi(6,k,i,j))*.5
        call freac(pol,cona,conm,sni)
c . . . Дополнительный джоуль от флуктуаций Эл. поля kiss=1
        kiss=0
        call djoulp (dj,vimm,vimd,vimr,viam,viad,viar,
     *               vin,vjn,vrn,ano,and,antr,sni,kiss)
c . . . Перекос в Джоуле для зимы в северном полушарии
        i_pol=8
        if (i.le.i_pol) then
          qdj(k)=dj  
	!	qdj(k)=2.*dj        ! НАГРЕВ ОДИНАКОВ
        else
          qdj(k)=dj  
	!	qdj(k)=2.*dj
        end if
        q(k)=e_dis*di+qi-co-cik53+che
    1 continue
      q(nh)=0.
      return
      end


C
      function ano ik(con no,ano,temp)
      real k17
c
      k17=6.5e-11
c
      a10=12
      rez=ano*k17
      w=a10/(a10+rez)
      anoik=3.73e-13*w*con no *rez*exp(-2700./temp)
      return
      end
c
      function cold(ano,and,antr,tem,anco2,radsk,gk,ctdk)
      real ltau,ksy
      data am1,am2,am3,amco2/53.12e-24,46.51e-24,26.56e-24,
     *     73.04e-24/,
     *     bk/1.38e-16/,re/6.371e8/
c*
      a=exp(-228/tem)
      f=1.67e-18*a*antr
      co o=f/(1.+0.6*a+0.2*exp(-325/tem))
c*
      d co2 o2=1.31e-14*exp(-41./(tem**0.33))*tem
      d co2 n2=dco2o2/3.3
      d co2 o =1.5e-11 *exp(-800./tem)
      sigma=dco2o2*ano+dco2n2*and+dco2o*antr
      w=1.2/(1.2+sigma)
      hco2=bk*tem/(amco2*gk)
      tau=hco2*anco2*6.4e-15
      l tau=e(tau)+e(tau/2.)
      ksy=0.5*w*l tau
      f ksy=ksy/(1.-w+ksy)
      co2ik=1.33e-13*2.*exp(-960./tem)*anco2*
     *       sigma*fksy
      cold=co2ik+co o
      return
      end
c
c
      subroutine djoul p(dj,vii,vij,vri,via i,via j,via r,
     *                   vi,vj,vr,an1,an2,an3,sni,kiss)
      dimension sni(6),ami(2),difv(2),am(3),an(3)
      data am/53.12,46.51,26.56/,
     *     ami/51.47,26.56/
      an(1)=an1
      an(2)=an2
      an(3)=an3
c
c . . . Увеличение действия электрического поля в шапке
c                           флуктуациями
      IF(kiss.eq.1) THEN
       const=1.5e4
       const=0.
       if(abs(vii).lt.const) then
          sig=sign(1.,vii)
          viif=vii+sig*const
       end if
       if(abs(vij).lt.const) then
          sig=sign(1.,vij)
          vijf=vij+sig*const
       end if
       if(abs(viai).lt.const) then
          sig=sign(1.,viai)
          viaif=viai+sig*const
       end if
       if(abs(viaj).lt.const) then
          sig=sign(1.,viaj)
          viajf=viaj+sig*const
       end if
       difv(1)=(viif-vi)**2+(vijf-vj)**2 +(vri-vr)**2
       difv(2)=(via if-vi)**2+(via jf-vj)**2 +(via r-vr)**2
      ELSE
       difv(1)=(vii-vi)**2+(vij-vj)**2 +(vri-vr)**2
       difv(2)=(via i-vi)**2+(via j-vj)**2 +(via r-vr)**2
      END IF
      dj=0.
      pm=0.
c*
      do 1 k=1,3
        pm=pm+sni(k)/((am(k)+ami(1))**2)
     *      * an(k)*am(k)
    1 continue
      pm=pm*ami(1)**2
      dj=pm*difv(1)*1.0e-24
c*
      pm=0.
      do 2 k=1,3
        pm=pm+sni(k+3)/((am(k)+ami(2))**2)
     *     * an(k)*am(k)
    2 continue
      pm=pm*ami(2)**2
      dj=dj+pm*difv(2)*1.0e-24
c*
      if(dj.lt.0.) go to 401
      return
  401 print 402,dj
  402 format(1x,'incorrect Joul heating',' dj=',e10.3)
  403 stop
      end
c
      subroutine freac(ant,cona,conm,sni)
      dimension sni(6)
      sni(1)=1.16e-9*conm
      sni(2)=1.50e-9*conm
      sni(3)=0.76e-9*conm
      sni(4)=1.21e-9*cona
      sni(5)=1.08e-9*cona
      sni(6)=1.86e-9*(ant/1000.)**0.37*cona
      return
      end
C
      function e(tau)
      real tau
      if(.not.((tau.ge.0.).and.(tau.le.10.)))
     *   go to 1
         e=0.5-4.76e-2*tau+2.67e-3*tau**2
         return
    1 continue
      if(.not.((tau.gt.10.).and.(tau.le.20.)))
     *   go to 2
         e=0.351-0.6e-2*tau
         return
    2 continue
      if(.not.(tau.gt.20.)) go to 3
        e=4.61/tau
        return
    3 continue
      print 100,tau
  100 format (2x,'function E: ПРОВЕРЬ tau','tau',e10.3)
      stop
      end
c
c          heatting term due to chem
c          reaction
c                                 o+o+m
      function chem(ano,and,antr,tem)
      a=300./tem
      f=8.16*4.7*1.e-22
      sum=ano+and+antr
c
      b=1.66e-23*exp(505./tem)
      chem= antr*1.e-23*sum*(f*a*antr+b*ano)
      if(chem.lt.1.e-15) chem=0.
c
      return
      end
c
