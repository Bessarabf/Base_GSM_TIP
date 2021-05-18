! subroutine terpot_bas, bongl, r_msis, co2con, nachc, nts, conn, progjn,
!      pgl3d, plotn, plots, boskli, bospgl, rezam, tnalt, turbk, cntn90, sumro
!      lowgln_bas, botcalc_L, noznew, intpa, bonPGL1, barsos, connot, timol, 
!      r_dis
!  function g11t31

c   terpot_bas - bas variant GSM TIP 2018-2019

c   ver.    
c   ver.    
c   ver.    
c   version 25.05.12 add to intrface KPA & NT for massive pril

c - conno1
c - sumro
      subroutine terpot_bas(day,god,dayt,godt,uts,tau,dts,solet,sole,
     *       solu,nsu,nse,kpars,rads,nh,gkoor,its,ddolgs,dtets,fa,fs,
     *       ap,pkp,dst,ae,al,au,bmpz,bmpy,mass,delta,pgl,pgi,ids,ins
     *      ,isp,vir,vid,vim,verno,parj,potef,ntr,nl2,pril,KPA,NT)
!   constants of GSM TIP   
      USE mo_bas_gsm
!   !!!!!!!!!!!!!!!!!!!!
      integer day,god,dayt,godt,verno
      dimension
     *   sole(nse),solu(nsu),solet(nse),parj(nh,its,ids)
     *   ,rads(nh),gkoor(2,its,ids),mass(30)
     *   ,pgi(ins,nh,its,ids),pgl(kpars,nh,its,ids)
     *   ,vir(nh,its,ids),vid(nh,its,ids),vim(nh,its,ids)
     *   ,potef(ntr,ids,nl2),pril(*)
     *   ,ros00(its0,ids0,nh0)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       allocatable an1(:,:,:),an2(:,:,:),an3(:,:,:) 
     *   ,an6(:,:,:)
     *   ,an11(:,:,:),an21(:,:,:),an31(:,:,:)
     *   ,an61(:,:,:),vr(:,:,:),vi(:,:,:)
     *   ,vj(:,:,:),vi1(:,:,:),vj1(:,:,:)
     *   ,g(:),rp(:),ctd(:),anco2(:,:,:)
     *   ,ron(:,:,:),ros(:,:,:),ros0(:,:,:)
    
       allocate (an1(its,ids,nh),an2(its,ids,nh),an3(its,ids,nh) 
     *   ,an6(its,ids,nh)
     *   ,an11(its,ids,nh),an21(its,ids,nh),an31(its,ids,nh)
     *   ,an61(its,ids,nh),vr(its,ids,nh),vi(its,ids,nh)
     *   ,vj(its,ids,nh),vi1(its,ids,nh),vj1(its,ids,nh)
     *   ,g(NH),rp(NH),ctd(NH),anco2(its,ids,nh)
     *   ,ros0(ITS,IDS,NH),ros(ITS,IDS,NH),ron(ITS,IDS,NH))

      data key/1/

      n4=nh

      do 21 k=1,nh
       g(k)=g0/(1.+rads(k)/re)**2
       if(k.ne.nh) then
         rp(k)=rads(k+1)-rads(k)
       else
         rp(nh)=rp(nh-1)*1.1
       end if
 21   continue
c
      call lowgln_bas(pgl,rads,kpars,nh,its,ids,day
     *         ,ap,fa,fs,gkoor,dtets,ddolgs,uts,mass(18),pril,KPA,NT)

       call pgl3d(pgl,kpars,nh,its,ids,an1,an2,an3,an6,vr,vi,vj)
!      rate dissociation massiv
       call r_dis(qdis,an1,an6,gkoor,g,rads,solu,nsu,delta,
     *            nh,its,ids,uts)
!      recommend time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dLong=(re+rads(nh-1))*sin(pi*dtets/180.)
      dLong=2*pi*dLong*ddolgs/360.
      Stau=dts+dts ! at the first time
      Vmax=abs(vj(its-1,1,nh))
      do j=2,ids
        if(Vmax.lt.abs(vj(its-1,j,nh))) Vmax=abs(vj(its-1,j,nh))
      end do
      if (Vmax.gt.0.1)Stau=dLong/Vmax
      if(Stau.lt.dts) then
         print*, 'GSMTIP: termos-recommended Stau=',Stau,Vmax,dLong
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call timol(pgl,kpars,nh,its,ids,dts,vim,vir,vid)
c
      call nachc(an1,an2,an3,an6,vi,vj,an11,an21,an31,an61
     *          ,vi1,vj1,nh,its,ids)
!!!      First step on time
      if(key.eq.1) then
        if((mass(5).eq.0).or.(mass(5).eq.2)) then
          call sumro(an1,an2,an3,ros,nh,its,ids)
        else
          call plots(ros,an1,an2,an3,an6,rp,g,nh,its,ids)
        end if
!       ROS0 ON THE FIRST STEP=ROS
        if((its.gt.its0).or.(ids.gt.ids0).or.(nh.gt.nh0)) then
          print*,' array ROS0 boundary exceed'
        stop
        end if
        ros0(1:its,1:ids,1:nh)=ros(1:its,1:ids,1:nh)
        key=0
      else
        ros0(1:its,1:ids,1:nh)=ros00(1:its,1:ids,1:nh)
        if((mass(5).eq.0).or.(mass(5).eq.2)) then
          call sumro(an1,an2,an3,ros,nh,its,ids)
        else
          call plots(ros,an1,an2,an3,an6,rp,g,nh,its,ids)
        end if
      end if
      call turbk   (ctd,rads,nh)
      ki=0
c
c     . . . NO-block
      if(mass(21).eq.0) then  ! apprpoximation
        call connot(pgl,rads,kpars,nh,its,ids)
      else
        call nonew(pgl,pgi,gkoor,ctd,rads,rp,g,
     *              kpars,ins,nh,its,ids,delta,uts,dts,mass)
      end if
c     . . . Расчет термосферы по MSIS
      if(mass(4).eq.0.or.mass(5).eq.0) then
         call r_msis(an11,an21,an31,an61,ron
     *                 ,rads,gkoor,mass
     *                 ,day,uts,nh,its,ids,fa,fs,ap)
      end if
c     . . . Температура рассчитывается
      if(mass(4).ne.0) then
         call co2con(anco2,an1,an2,an3,an6,ctd,
     *               rads,rp,g,nh,its,ids)
         call heatpo_bas(pgl,pgi,parj,solet,solu,nsu,nse,
     *               kpars,rads,g,nh,gkoor,its,ddolgs,dtets,
     *               mass,delta,day,uts,tau,dts,ctd,vim,vid,
     *               vir,ids,ins,an1,an2,an3,anco2,an6,an61,ros,
     *               vi,vj,vr)


         call nts (an61,nh,its,ids,nh,its-3)

      end if

      if(mass(5).ne.0) then
         call sdizkn_bas(an1,an2,an3,an11,an21,an31,
     *               an61,vr,vi,vj,ros,rp,rads,g,nh,its,ids,
     *               dts,ctd,roS,solu,gkoor,delta,nsu,
     *               dtets,uts,ddolgs)
         call nts (An11,nh,its,ids,nh,its-3) ! ,1
         call nts (An31,nh,its,ids,nh,its-3) ! ,1)  
      end if
c
      if(mass(5).eq.1) then
         call plots(ron,an11,an21,an31,an61,rp,g,nh,its,ids)
	  ELSE
	   call sumro(an11,an21,an31,roN,nh,its,ids)
      end if
c     . . . V = 0
    9 if(mass(6).ne.0) then
         call veter_ham(vi1,vj1,vi,vj,vr,vim,vid,an1,an2,an3,an61,ron,
     *              pgl,pgi,rp,rads,ctd,nh,its,ids,kpars,ins,dts)
         call nts (vj1,nh,its,ids,nh,its-3)
         call nts (vi1,nh,its,ids,nh,its-3)
         call bonvec1(vi1,vj1,nh,its,ids)
      else
         vi1=0.
         vj1=0.
         vi=0.
         vj=0.

      end if
c     . . . Vr = 0 ?
      if(mass(7).ne.0) then
	
        call vrprim_s(vr,vi,vj,vi1,vj1,ros0,ros,ron,rp,rads,
     *                nh,its,ids,dts,mass(10))
        call nts (vr,nh,its,ids,nh,1) ! its-3)!  1)
        call boskli(vr,nh,its,ids)
        call tnalt(vr,nh,its,ids,mass(10))
	
      else
        vr=0.
      end if
!    
c     ROS0 
      ros00(1:its,1:ids,1:nh)=ros(1:its,1:ids,1:nh)
!
      call nachc(an11,an21,an31,an61,vi1,vj1,
     *           an1,an2,an3,an6,vi,vj,nh,its,ids)

      call rezam(pgl,an11,an21,an31,an61,vi1,vj1,
     *           vr,kpars,nh,its,ids)
c          File: labt.dan writing heat sourse
      if(mass(19).ne.0) then
        open(8,file='labt.dan',form='Unformatted')
        rewind8
        call heapot(pgl,pgi,an11,an21,an31,an61,anco2,vr,
     *              vi1,vj1,vim,vid,vir,ctd,solu,
     *              nsu,rads,g,gkoor,kpars,ins,
     *              nh,its,ids,delta,uts,mass(19))
        close(8)
      end if
      print 137,mass(10)
  137 format(' thermospher calculate until',i4,
     *       ' altitude point')
	
      deallocate (an1,an2,an3
     *   ,an6,an11,an21,an31
     *   ,an61,vr,vi
     *   ,vj,vi1,vj1
     *   ,g,rp,ctd,anco2
     *   ,ron,ros,ros0)

      return
      end
!------------------------------------------------------------------------------
      subroutine bongl(an1,n,n1,n2)
      dimension an1(n1,n2,n)
c 900 format(' bongl  :',3i7)
c      print 900,n,n1,n2
      np=n1-1
      n6=n2/2
      do1k=1,n
      do2j=1,n6
      j1=j+n6
      an1(1,j,k)=(an1(2,j,k)+an1(2,j1,k))/2.
      an1(1,j1,k)=an1(1,j,k)
      an1(n1,j,k)=(an1(np,j,k)+an1(np,j1,k))/2.
      an1(n1,j1,k)=an1(n1,j,k)
   2  continue
  1   continue
      return
      end
!------------------------------------------------------------------------------
c     . . . Расчет термосферы по MSIS
      subroutine r_msis(an11,an21,an31,an61,ron
     *                 ,rads,gkoor,mass
     *                 ,day,uts,nh,its,ids,fa,fs,ap)
      dimension an11(its,ids,nh),an21(its,ids,nh),an31(its,ids,nh),
     *          an61(its,ids,nh),ron(its,ids,nh),
     *          rads(nh),gkoor(2,its,ids)
      dimension cns1(30),cns2(30),cns3(30),tns(30),mass(30)
      integer day,god,dayt,godt
      data pi,om/3.1415926,7.27e-5/,
     *     am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
      do 88 i=1,its
        do 88 j = 1,ids
           fig=gkoor(1,i,j)
           fig=90.-fig
           dgeo=gkoor(2,i,j)/180.*pi
           dol=gkoor(2,i,j)
           tau1=uts+dgeo/om
           tau1=tau1/3600.
           tau2=uts
           if(tau2.gt.86400.) tau2=tau2-86400.
           td=day
           if(mass(17).eq.86) then
             call cntn90(rads,fig,dol,td,ap,fa,fs,
     *                  tau1,tau2,cns2,cns1,cns3,tns,nh)
             do k=1,nh
                an61(i,j,k)=tns(k)
                if(mass(5).eq.0) then
                  an11(i,j,k)=cns1(k)
                  an21(i,j,k)=cns2(k)
                  an31(i,j,k)=cns3(k)
c                 an31(i,j,k)=cns3(k)/1.4
                  ron(i,j,k)=an11(i,j,k)*am1+an21(i,j,k)*am2
     *                       +an31(i,j,k)*am3
                end if
            end do
           else
              print 884,mass(17)
              stop
           end if
   88    continue
  884    format('  terpot: mass(17)=',i4,'  Program stop')
         return
         end
!------------------------------------------------------------------------------
      subroutine co2con(an co2,an1,an2,an3,an6,ctd,
     *                   rads,rp,g,n,n1,n2)
      dimension an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n),ctd(n)
     *         ,an6(n1,n2,n),anco2(n1,n2,n),rads(n),g(n),rp(n)
      data am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/,
     *     amco2/73.04e-24/,bk/1.38e-16/
c
      do 11 k=1,n
       z=rads(k)*1.e-5
       if(.not.(z.gt.100.)) go to 12
        k1=k
        goto 13
   12  continue
   11 continue
   13 do 1 i=1,n1
       do 2 j=1,n2
        do 3 k=1,k1
         sum con=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
c
         ams=(am1*an1(i,j,k)+am2*an2(i,j,k)+am3*
     *        an3(i,j,k))/sum con
          z km=rads(k)*1.e-5
          h cm=bk*an6(i,j,k)/(ams*g(k))
          hkm=hcm*1.e-5
          anco2(i,j,k)=1.e11*an6(i,j,1)/an6(i,j,k)*
     *              exp((1.+sqrt(1.+0.8e-12*hcm*hcm))*
     *              (80.-zkm)/(2.*hkm))
    3    continue
    2  continue
    1 continue
      call barsos(anco2,an6,rp,g,amco2,n,n1,n2,k1)
      return
      end
!------------------------------------------------------------------------------
      subroutine conn(pgl,kpars,nh,its,ids)
       dimension pgl(kpars,nh,its,ids)
          do 1 k = 1 , nh
          do 1 i = 1 , its
          do 1 j = 1 , ids
       pgl(5,k,i,j)=0.
  1     continue
       return
       end
!------------------------------------------------------------------------------
      subroutine nachc(an1,an2,an3,an6,vi,vj,
     *                 an11,an21,an31,an61,vi1,vj1,n,n1,n2)
c
c     reevaluate array after the time-step
c
      dimension an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n),
     *          an6(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n),
     *          an11(n1,n2,n),an21(n1,n2,n),an31(n1,n2,n),
     *          an61(n1,n2,n),vi1(n1,n2,n),vj1(n1,n2,n)
         an11=an1
         an21=an2
         an31=an3
         an61=an6
         vi1=vi
         vj1=vj
      return
      end
!------------------------------------------------------------------------------
       subroutine nts(an1,nn,n1,n2,n,istep)
!      zonal smoothing
!      beta - smoothing parameter
!      istep - latitude step
      dimension an1(n1,n2,nn)
      allocatable pm(:),pmm(:),ap(:),bp(:),cp(:),fp(:)
      allocate(pm(n2),pmm(n2),ap(n2),bp(n2),cp(n2),fp(n2))
      beta=0.15
      ni=n1-1
       do  k=2,n
          do i=2,n1-1,istep
	       do j=1,n2
c***********
             dd=0.
             ap(j)=beta
             bp(j)=1.+2.*beta
             cp(j)=beta
!           fp(j)=-alog(an1(i,j,k))
             fp(j)=-an1(i,j,k)
           end do
           call progjn(ap,bp,cp,fp,n2,np,pm,pmm)
           do j=1,n2
              an1(i,j,k)=(pm(j))
           end do
          end do
        end do
        deallocate(pm,pmm,ap,bp,cp,fp)
        return
        end
!------------------------------------------------------------------------------
	subroutine progjn(aa,bb,cc,ff,n2,np,pm,pmm)
      dimension aa(n2),bb(n2),cc(n2),ff(n2),pm(n2),pmm(n2)
      allocatable alf(:),bet(:),u(:),vp(:)
     *            ,e(:),d(:),gg(:)
      allocate (alf(n2),bet(n2),u(n2),vp(n2)
     *         ,e(n2),d(n2),gg(n2))
       n3=n2-1
       np=n2+1
       e(1)=cc(1)/bb(1)
       d(1)=aa(1)/bb(1)
       gg(1)=-ff(1)/bb(1)
       do 1 i=2,n2
         ab=bb(i)-aa(i)*e(i-1)
         e(i)=cc(i)/ab
c
         d(i)=aa(i)*d(i-1)/ab
         gg(i)=(aa(i)*gg(i-1)-ff(i))/ab
   1   continue
       ac=bb(n2)-aa(n2)*(e(n3)+d(n3))
       alf(n2)=cc(n2)/ac
       bet(n2)=(aa(n2)*gg(n3)-ff(n2))/ac
       do 2 i=1,n3
         j=n2-i
         alf(j)=e(j)*alf(j+1)+d(j)*alf(n2)
         bet(j)=e(j)*bet(j+1)+d(j)*bet(n2)+gg(j)
   2   continue
       pm(1)=(-ff(1)+aa(1)*bet(n2)+cc(1)*bet(2))
     *      /(bb(1)-aa(1)*alf(n2)-cc(1)*alf(2))
       pmm(n2)=alf(n2)*pm(1)+bet(n2)
       do 3 i=2,n2
        pm(i)=bet(i)+pm(1)*alf(i)
  3    continue
       do 4 i=1,n3
         j=n2-i
         pmm(j)=e(j)*pmm(j+1)+d(j)*pmm(n2)+gg(j)
 4     continue
       deallocate (alf,bet,u,vp
     *            ,e,d,gg)
       return
       end
!------------------------------------------------------------------------------
      subroutine pgl3d(pgl,kpars,nh,its,ids,
     *                  an1,an2,an3,an6,vr,vi,vj)
      dimension pgl(kpars,nh,its,ids),an1(its,ids,nh)
     *       ,an2(its,ids,nh),an3(its,ids,nh),an6(its,ids,nh)
     *       ,vr(its,ids,nh),vi(its,ids,nh),vj(its,ids,nh)
        do 1 j=1,ids
         do 2 i = 1 , its
           do 3 k = 1 , nh
           an1(i,j,k)=pgl(1,k,i,j)
           an2(i,j,k)=pgl(2,k,i,j)
           an3(i,j,k)=pgl(3,k,i,j)
           an6(i,j,k)=pgl(7,k,i,j)
           vr(i,j,k)=pgl(10,k,i,j)
           vi(i,j,k)=pgl(11,k,i,j)
           vj(i,j,k)=pgl(12,k,i,j)
  3       continue
   2     continue
    1   continue
      return
      end
!------------------------------------------------------------------------------
c      28.02.96 новая плотность
      subroutine plotn(ro,an1,an2,an3,an6,r,
     *                 g,n,n1,n2)
      dimension ro(n1,n2,n),an1(n1,n2,n),an2(n1,n2,n),g(n),
     *          an3(n1,n2,n),an6(n1,n2,n),r(n)
      data am1,am2,am3/53.12e-24,46.56e-24,26.56e-24/
     
      data rg/8.31442e07/
      n11=n1-1
      do 1 i=2,n11
       do 2 j=1,n2
         ro(i,j,1)=am1*an1(i,j,1)+am2*an2(i,j,1)+
     *             am3*an3(i,j,1)
         do 3 k=2,n
          dx=r(k)-r(k-1)
          ano=an1(i,j,k-1)
          and=an2(i,j,k-1)
          antr=an3(i,j,k-1)
          amcn=32.*ano+28.*and+16.*antr
          amcn=amcn/(ano+and+antr)
          ano=an1(i,j,k)
          and=an2(i,j,k)
          antr=an3(i,j,k)
          amcv=32.*ano+28.*and+16.*antr
          amcv=amcv/(ano+and+antr)
          tn=an6(i,j,k-1)
          tv=an6(i,j,k)
         rt12=rg*(1./amcv+1./amcn)/2.*(tv+tn)/2.
         aa=(g(k)+g(k-1))/2.+
     *      (tn+tv)/2.*rg*(1./amcv-1./amcn)/dx+
     *      rg*(1./amcv+1./amcn)/2.*(tv-tn)/dx
         ro(i,j,k)=ro(i,j,k-1)*(rt12/dx-aa/2.)/
     *               (rt12/dx+aa/2.)
    3    continue
    2  continue
    1 continue
      call bongl(ro,n,n1,n2)
      return
      end
!------------------------------------------------------------------------------
      subroutine plots(ro,an1,an2,an3,an6,rp,
     *                g,n,n1,n2)
c      7.06.98 новая плотность, уточненная
      dimension ro(n1,n2,n),an1(n1,n2,n),an2(n1,n2,n),g(n),
     *          an3(n1,n2,n),an6(n1,n2,n),rp(n)
      data am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/,
     *     bk/1.38e-16/
      do 1 i=1,n1
       do 2 j=1,n2
         ro(i,j,1)=am1*an1(i,j,1)+am2*an2(i,j,1)+
     *             am3*an3(i,j,1)
	 
         do 3 k=2,n
	
          ano=an1(i,j,k-1)
          and=an2(i,j,k-1)
          antr=an3(i,j,k-1)
          amcn=am1*ano+am2*and+am3*antr
	
          amcn=amcn/(ano+and+antr)
          ano=an1(i,j,k)
          and=an2(i,j,k)
          antr=an3(i,j,k)
          amcv=am1*ano+am2*and+am3*antr
          amcv=amcv/(ano+and+antr)
c
          tn=an6(i,j,k-1)
          tv=an6(i,j,k)
c . . . средняя шкала высот (обратная)
          oh1=g(k-1)*amcn/bk/tn
          oh2=g(k)*amcv/bk/tv
          alf=(oh1+oh2)*rp(k-1)*0.5
          s=alog(ro(i,j,k-1)*tn*amcv/(tv*amcn))-alf
          ro(i,j,k)=exp(s)
    3    continue
    2  continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine boskli(an1,nh,its,ids)
C . . . Расчет в полюсной точке. Посвящается Клименко В.И.
C   . . . ids - четное
      dimension an1(its,ids,nh)
      i2=its-1
      do 1 k=1,nh
       s sp=0.
       s np=0.
       do 2 j=1,ids
        s np=snp+an1(2,j,k)
        s sp=s sp+an1(i2,j,k)
    2  continue
       u np=s np/ids
       u sp=s sp/ids
       do 3 j=1,ids
        an1(1,j,k)=u np
        an1(its,j,k)=u sp
    3  continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
	subroutine bospgl(pgl,kpars,nh,its,ids,np)
c . . . Расчет в полюсной точке. Посвящается Клименко В.И.
C   . . . ids - четное
      dimension pgl(kpars,nh,its,ids)
      i2=its-1
      do 1 k=1,nh
c     ssp - sum s.pole
c     snp - sum n.pole
       s sp=0.
       s np=0.
       do 2 j=1,ids
        s np=snp+pgl(np,k,2,j)
        s sp=s sp+pgl(np,k,i2,j)
    2  continue
       u np=s np/ids
       u sp=s sp/ids
       do 3 j=1,ids
        pgl(np,k,1,j)=u np
        pgl(np,k,its,j)=u sp
    3  continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine rezam(pgl,an1,an2,an3,an6,vi,vj,vr,kpars,nh,its,ids)
      dimension pgl(kpars,nh,its,ids),an1(its,ids,nh),an2(its,ids,nh)
     *,an3(its,ids,nh),an6(its,ids,nh),vi(its,ids,nh),vj(its,ids,nh)
     *,vr(its,ids,nh)
        do 1 k = 1 , nh
        do 1 i = 1 , its
        do 1 j = 1 , ids
        pgl(1,k,i,j)=an1(i,j,k)
        pgl(2,k,i,j)=an2(i,j,k)
        pgl(3,k,i,j)=an3(i,j,k)
        pgl(7,k,i,j)=an6(i,j,k)
        pgl(11,k,i,j)=vi(i,j,k)
        pgl(12,k,i,j)=vj(i,j,k)
        pgl(10,k,i,j)=vr(i,j,k)
  1     continue
        return
        end
!------------------------------------------------------------------------------
      subroutine tnalt(an6,n,n1,n2,n0)
      dimension an6(n1,n2,n)
      nm=n-1
      do 1 k=n0,nm
       do 2 j=1,n2
        do 3 i=1,n1
         an6(i,j,k+1)=an6(i,j,k)
  3     continue
  2    continue
  1   continue
      return
      end
!------------------------------------------------------------------------------
      subroutine turbk(ctd,rads,n)
      dimension ctd(n),rads(n)
c      data hm/103./,c1,c0/1.e 6,5.e 5/,
c      data hm/ 94./,c1,c0/1.e 6,5.e 5/,
c    *     s1,s2,s3/0.05,0.05,0.07/
c      . . . Вариант с увеличенным Кт и медленным спадом
c  !     data hm/ 94./,c1,c0/5.e 6,5.e 6/,
c      data hm/ 94./,c1,c0/1.e 6,5.e 5/,
!      data hm/ 90./,c1,c0/5.e 6,1.e 6/,!9.03.11
!     *     s1,s2,s3/0.01,0.05,0.01/
	data hm/98./,c1,c0/5.0e 6,3.0e 6/
     *     s1,s2,s3/0.01,0.005,0.01/
!	data hm/103./,c1,c0/5.0e 6,3.0e 6/	 !!! 5/02/13
!        data hm/108./,c1,c0/8.0e 6,4.0e 6/	 !!! 25/07/18
!     *     s1,s2,s3/0.01,0.005,0.01/
      do k=1,n
       ah=(rads(k))*1.e-5
       y=hm-ah
       y2=y*y
       if (y2.gt.1000.) then
         ctd(k)=0.
       else
         if(ah.lt.hm) then
           ctd1=(c1-c0)*exp(-s2*y2)
           ctd(k)=ctd1+c0*exp(-s3*y)
         else
           ctd(k)=c1*exp(-s1*y2)
         end if
       end if
      end do
      return
      end
!------------------------------------------------------------------------------
c        Пoдпpoгpaммы для  MSIS - 86 !!!
c    Головная названа CNTN, чтобы вызывалась одинаково
c        с другими версиями ( MSIS-77, MSIS-83 ).
c
c   Автоответчик версии:
c     function nmsis
c     nmsis=86
c     return
c     end
c
c
c  П/п paccчитывaeт пapaмeтpы нeйтpaльнoй aтмocфepы пo мoдeли MCИC-90
!------------------------------------------------------------------------------
      subroutine cntn90(rads,fig,dolg,td,ap,fa,fs,tlttau,utsec,
     *                 cn1,cn2,cn4,tn,nh)
c   cn4,cn2,cn1,dh,tn - мaccивы концентраций O,O2,N2,H и тeмпepaтyp
c    в yзлax на высотах rads в точке fig,dolg в момент td,uttau
      dimension rads(nh),cn1(nh),cn2(nh),cn4(nh),tn(nh)
      dimension t(2),apm(7),d(9) ! MSIS2000
!     dimension d(8)             ! MSIS86
      do 3 i=1,7
       apm(i)=ap   ! ap-index
 3    continue

      alat=fig*1.7453292e-2
      alon=dolg*1.7453292e-2
      iyd=80*1000+td
     
      do 2 ih=1,nh
!           MSIS-86
!            call gtd6(iyd,utsec,rads(ih)/1.e5,fig,dolg,tlttau,fs,fa,apm
!     *      ,48,d,t)
!           MSIS 2000
         call gtd7(iyd,utsec,rads(ih)/1.e5,fig,dolg,tlttau,fs,fa,apm
     *            ,48,d,t)
         
         tn(ih)=t(2)
         cn1(ih)=d(3)
         cn2(ih)=d(4)
         cn4(ih)=d(2)
 2    continue

      return
      end
!------------------------------------------------------------------------------
      subroutine sumro(an1,an2,an3,ro,n,n1,n2)
      dimension an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n)
     *         ,ro(n1,n2,n)
      data am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
      do 1 i=1,n1
       do 2 j=1,n2
        do 3 k=1,n
         ro(i,j,k)=an1(i,j,k)*am1
     *            +an2(i,j,k)*am2
     *            +an3(i,j,k)*am3
    3  continue
    2  continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
c . . . ver. 2012 - P_riliv 
      subroutine lowgln_bas(pgl,rads,kpars,nh,its,ids,day
     *           ,ap,fa,fs,gkoor,dtets,ddolgs,uts,musl,pril,KPA,NT)
c
c     . . . параметры на 80 км
c
!      USE mo_ham_gsm 
      
      dimension pgl(kpars,nh,its,ids),rads(nh)
     *         ,gkoor(2,its,ids)
     *         ,apm(7),tm(2)
!      dimension dm(8) ! MSIS 90
      dimension dm(9) ! MSIS2000


      dimension pril(*)
      integer day
      data om,pi/7.27e-5,3.14159/
     *    ,am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
     *     ,bk,re,gg/1.38e-16,6.371e8, 956.81976/
      i1=2
      i2=its-1
      nrm=kpars
      npt=7
       do 100 i=1,7
 100       apm(i)=ap
c     . . .  NO and N
      IF(musl.NE.3) THEN
      do i=1,its
       do j=1,ids
        pgl(4,1,i,j)=1.e+06
        pgl(5,1,i,j)=5.e+04
       end do
      end do
      END IF
      IF(musl.eq.0) THEN
      do 1 i =  1,its
       do 1 j = 1,ids
c      . . . constant value on lower boundary
c      . . .  m usl=0 - yes
c       fig=gkoor(1,5,1)
c       fig=90.-fig
c       dgeo=gkoor(2,5,1)/180.*pi
c       tau=12
c       td=day
C . . .  F10.7=70
        pgl(1,1,i,j)=7.4e+13
        pgl(2,1,i,j)=3.0e+14
C . . .  F10.7=180
c        pgl(1,3,i,j)=2.1e+13
c       pgl(2,1,i,j)=3.2e+14
c
        pgl(3,1,i,j)=2.4e+11
        pgl(7,1,i,j)=180.
cc        pgl(3,1,i,j)=2.4e+11
cc        pgl(7,1,i,j)=180.

   1   continue
       ELSE IF(m usl.eq.1) THEN
c    . . . lowboundary on MSIS-90
        do   i =  1,its
          do   j = 1,ids
            fig=gkoor(1,i,j)
            fig=90.-fig
            dgeo=gkoor(2,i,j)/180.*pi
            dol=gkoor(2,i,j)
            tau1=uts+dgeo/om
            tau1=tau1/3600.
c . . . Сдвиг фазы на 2 часа
            tau3=tau1-2.
            tau2=uts-7200.
            if(tau2.gt.86400)tau2=tau2-86400.
            td=day
c  12       alt=rads(1)/1.e5
            iyd=80*1000+td
            vis=rads(1)/1.e5
!            call gtd6(iyd,tau2,80.,fig,dol,tau3,fa,fs,
!     *                 apm,48,dm,tm)
            call gtd7(iyd,tau2,80.,fig,dol,tau3,fa,fs,
     *                 apm,48,dm,tm)
            pgl(7,1,i,j)=tm(2)
            pgl(1,1,i,j)=dm(4)
            pgl(2,1,i,j)=dm(3)
            pgl(3,1,i,j)=dm(2)
           end do
          end do
          ELSE IF(m usl.eq.2) THEN
               td=day
               call botcalc_L(pgl,nh,its,ids,kpars,uts,dtets,ddolgs,
     *                      rads,gkoor,pril,KPA,NT)
!         HAMMONIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	         
!          ELSE IF(m usl.eq.3) THEN
!            do i =  1,its
!              do  j = 1,ids 
!                pgl(7,1,i,j)=gsmHAM(1,i,j)
!
! sum concentration
!                aNall=dgsmHAM(1,i,j)/(0.21*am1+0.79*am2)
!
!                pgl(1,1,i,j)=0.21*aNall   !!! 05.03.19
!                pgl(2,1,i,j)=0.79*aNall
!
!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!
!                pgl(3,1,i,j)=dOgsm(1,i,j)
!                
!                pgl(4,1,i,j)=dNOgsm(1,i,j)
!                pgl(5,1,i,j)=dNgsm(1,i,j)
!!!!!!!!!!!!!!!!!!!!!! 
!               pgl(11,1,i,j)=UgsmHAM(1,i,j)
!                pgl(12,1,i,j)=VgsmHAM(1,i,j)
!
!              end do
!            end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
          ELSE
                print *,'GSMTIP: lowgln incorrect МАSS(18) in lowgln'
                stop
c               END IF
           END IF
      return
      end
!------------------------------------------------------------------------------
!         Считывание данных Miyahara (или других)
      subroutine botcalc_L(pgl,nh,its,ids,kpars,uts,dtets,ddolgs,
     *           rads,gkoor,pril,kpa,ntime)
     
     
      dimension pgl(kpars,nh,its,ids),gkoor(2,its,ids),rads(nh)
   
      dimension pril(kpa,its,ids,ntime),mkp(7)

      allocatable p(:,:,:),gir(:),tes(:)

      data mkp/1,2,3,7,10,11,12/
      allocate (p(kpars,its,ids),gir(ids),tes(its))

  900 format (a25)
  901 format (a1)
  902 format (F5.1,4F8.2)
  903 format (F5.1,3X,E9.2,8X,E9.2,5X,3F8.1)
  904 format (5E11.3)
!      расчет момента времени
!     данные записаны с интервалом ntime час
	nDT=24/ntime        ! time interval
!      l=(uts-dts)/3600./nDT+1   ! number point of massive PRIL
      l=(uts-0.1)/3600./nDT+1   ! number point of massive PRIL
      if(l.gt.ntime) then
		print*, 'massive pril exeeded in botcalc ',l,' uts=',uts
!	stop
		l=ntime
      end if
!  longitude
      do i=1,ids
        gir(i)=ddolgs*(i-1)
      enddo
!  latitude 
      do i=1,its
        tes(i)=dtets*(i-1)
      enddo
C  для концентраций                
      CN=2.96e14+7.95e13+8.5e10 ! summury density
      A1=7.95E13/CN                
      A2=2.96E14/CN                
      A3=8.5E10/CN                 
      AMP=48.12E-24                
C**********************************************************
      do j = 1,ids
        do i = 1,its
	     prr=pril(1,i,j,l) !
            p(1,i,j)=prr*a1/amp ! O2
           p(2,i,j)=prr*a2/amp ! N2
           p(3,i,j)=prr*a3/amp ! O
           p(7,i,j)=pril(2,i,j,l)        ! T
	     p(10,i,j)=pril(3,i,j,l)       !
           p(11,i,j)=pril(4,i,j,l)       !
           p(12,i,j)=pril(5,i,j,l)       !
        enddo
      enddo
      do ikp=1,7 ! interpolation
       kpar=mkp(ikp) 
       call intpa(tes,its,dtets,gir,ids,ddolgs,gkoor,p,  ! interpolation to geom 
     *            kpar,pgl,kpars,nh)
      enddo
	
      call bonPGL1(pgl,kpars,nh,its,ids)
      
      call noznew( gir,tes,kpars,nh,its,ids,pgl) ! vector to geomag coor

	print*,'botcalc: end of interpolation'
      deallocate (p,gir,tes)
      return
      end
!------------------------------------------------------------------------------
      subroutine noznew(gir,tes,kpars,nh,its,ids,pgl)
      dimension gir(ids),pgl(kpars,nh,its,ids),tes(its)
       i=1
        do j=1,its
          do l=1,ids
            t=tes(j)
            f=gir(l)
            r=g11t31(f,t)
            r=-r
            s=sin(r)
            c=cos(r)
            st=pgl(11,1,j,l)
            sf=pgl(12,1,j,l)
            p1=st*c-sf*s
            pgl(11,1,j,l)=p1
            p2=sf*c+st*s
            pgl(12,1,j,l)=p2
          enddo
        enddo
      return
      end
!------------------------------------------------------------------------------
      subroutine intpa(tes,its,dtets,gir,ids,ddolgs,gkoor,p,
     *            kpar,pgl,kpars,nh)
      dimension tes(its),gir(ids),gkoor(2,its,ids),p(kpars,its,ids),
     *          pgl(kpars,nh,its,ids)
  900 format (2e12.3)
      do j=1,ids
        do i=1,its
          tet =gkoor(1,i,j)
          dolg=gkoor(2,i,j)
c         print *,' i,tet,dolg',i,tet,dolg
          call find(its,tet, tes,in)
          call find(ids,dolg,gir,jn)
          tetin=tes(in)
          dx=(tet-tetin)/dtets
          dolgjn=gir(jn)
c         print *,' tetin,dolgjn',tetin,dolgjn
          dy=(dolg-dolgjn)/ddolgs
          jn1=jn+1
          if(jn.eq.ids)jn1=1
c         print *,' in,in+1,jn,jn+1',in,in+1,jn,jn1
          p1=p(kpar,in,jn)
          p2=p(kpar,in+1,jn)
          p3=p(kpar,in,jn1)
          p4=p(kpar,in+1,jn1)
          f1=p1+(p2-p1)*dx
          f2=p3+(p4-p3)*dx
          ff=f1+(f2-f1)*dy
c         print *,' p1,p2,p3,p4,ff',p1,p2,p3,p4,ff
          pgl(kpar,1,i,j)=ff
        enddo
      enddo
      return
      end
!------------------------------------------------------------------------------
      subroutine bonPGL1(pgl,kpars,nh,its,ids)
      dimension pgl(kpars,nh,its,ids)
     *         ,inp(4)        
	data inp/1,2,3,7/
      i2=its-1
      k=1

	do 3 i=1,3
          np=inp(i)
         
c      ssp - sum s.pole
c      snp - sum n.pole
        s np=0.
        s sp=0.
        do 1 j = 1 , ids
         snp=snp+pgl(np,k,2,j)
         ssp=ssp+pgl(np,k,i2,j)
   1    continue
c
        unp=snp/ids
        usp=ssp/ids
        do 4 j=1,ids
          pgl(np,k,1,j)=unp
          pgl(np,k,its,j)=usp
    4   continue
    
    3 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine barsos(an,an6,rp,g,am,n,n1,n2,l)
      dimension an(n1,n2,n)
     *         ,an6(n1,n2,n),rp(n),g(n)
      data bk/1.38e-16/,ves/0.5/
      f2=bk/am
      do 1 i=1,n1
       do 2 j=1,n2
        do 3 k=l,n
          tn=an6(i,j,k-1)
          tv=an6(i,j,k)
          h1=f2*tn/g(k-1)
          h2=f2*tv/g(k)
          alf=(ves/h1+(1.-ves)/h2)*rp(k-1)
          ss=alog(an(i,j,k-1)*tn/tv)-alf
          an(i,j,k)=exp(ss)
    3   continue
    2  continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine bonvec1(vi,vj,n,n1,n2)
      dimension
     *         vi(n1,n2,n),vj(n1,n2,n)
      data pi/3.1415926/
      np=n1-1
      dfi=pi*2.0/n2
      dtet=pi/np
      n6=n2/2
      cosin=cos(dtet)
      do 1 k=1,n
       do 2 j=1,n6
        j1=j+n6
        fi=(j-1)*dfi
        fi180=fi+pi
c    . . . north pole
        as=sin(cosin*fi)
        ac=cos(cosin*fi)
        asin pi=sin(cosin*fi180)
        acos pi=cos(cosin*fi180)
        a1=vj(2,j,k)*ac+
     *     vi(2,j,k)*as
        a2=vj(2,j1,k)*acos pi+
     *     vi(2,j1,k)*asin pi
        a=0.5*(a1+a2)
        b1=vj(2,j,k)*as-vi(2,j,k)*ac
        b2=vj(2,j1,k)*asin pi-
     *     vi(2,j1,k)*acos pi
        b=(b1+b2)*0.5
        vi(1,j,k)=a*sin(fi)-b*cos(fi)
        vi(1,j1,k)=-vi(1,j,k)
        vj(1,j,k)=a*cos(fi)+b*sin(fi)
        vj(1,j1,k)=-vj(1,j,k)
c    . . . south pole
        a1=vj(np,j,k)*ac-
     *     vi(np,j,k)*as
        a2=vj(np,j1,k)*acos pi-
     *     vi(np,j1,k)*asin pi
        a=0.5*(a1+a2)
        b1=-vj(np,j,k)*as-vi(np,j,k)*ac
        b2=-vj(np,j1,k)*asin pi-
     *     vi(np,j1,k)*acos pi
        b=(b1+b2)*0.5
        vi(n1,j,k)=-a*sin(fi)-b*cos(fi)
        vi(n1,j1,k)=-vi(n1,j,k)
        vj(n1,j,k)=a*cos(fi)-b*sin(fi)
        vj(n1,j1,k)=-vj(n1,j,k)
   2   continue
   1  continue
      return
      end
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
      subroutine timol(pgl1,kpars,nh,its,ids,dts,vim,vir,vid)
      real vim(nh,its,ids),vid(nh,its,ids),vir(nh,its,ids),nu0,
     *     pgl1(kpars,nh,its,ids)
      data bk/1.38e-16/,nu0/0.9e-9/,ami/30./,ae/1.6e-24/
  900 format(' ',10g12.3)
      a0=2*7.69e-19/3/bk
      c0=ami*ae*nu0/6/bk
      igp=its-1
      do 3 j = 1 , ids
       do 1 ig=2,its-1
        do 2 i=1,nh
      if(pgl1(9,i,ig,j).le.0.) print22,pgl1(9,i,ig,j),i,ig,j
  22   format(' error , te=',e10.3,' i=',i4,' ig=',2i4)
          dvr=pgl1(10,i,ig,j)-vir(i,ig,j)
          dvt=pgl1(11,i,ig,j)-vim(i,ig,j)
          dvf=pgl1(12,i,ig,j)-vid(i,ig,j)
          dvr2=dvr*dvr
          dvt2=dvt*dvt
	

          dvf2=dvf*dvf
          a=a0*pgl1(6,i,ig,j)/(pgl1(9,i,ig,j)**1.5)
          sn=pgl1(1,i,ig,j)+pgl1(2,i,ig,j)+0.45*pgl1(3,i,ig,j)
          b=nu0*sn/2.
          c=c0*sn*(dvr2+dvt2+dvf2)
          r1=dts*(a*pgl1(9,i,ig,j)+b*pgl1(7,i,ig,j)+c)
          r2=1.+dts*(a+b)
          pgl1(8,i,ig,j)=(pgl1(8,i,ig,j)+r1)/r2
          if(pgl1(8,i,ig,j).le.pgl1(7,i,ig,j)) pgl1(8,i,ig,j)=
     *   pgl1(7,i,ig,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          pgl1(8,i,ig,j)=pgl1(7,i,ig,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	     
    2   continue
    1 continue
    3   continue
      return
      end
!------------------------------------------------------------------------------
      function g11t31(f,t)
      double precision a,b,c,e,ps,r,s,u,v
      data s/1.976573d-1/,c/9.802712d-1/
      data ps/1.74532925199432d-2/
      e=dble(t)*ps
      u=dble(f)*ps
      a=dcos(e)
      b=dcos(u)
      r=a*b*s
      a=dsin(e)
      v=c*a
      e=v+r
      a=dsin(u)
      r=s*a
      b=datan2(r,e)
      g11t31=sngl(b)
      return
      end
!------------------------------------------------------------------------------
c ver 20/03/2019
c 4d massive dissociation rates (1 - 1/cm3/s; 2 - erg/cm3/s )
      subroutine r_dis(qdis,ano2,tem,gkoor,g,rads,solu,nsu,del,
     *           nh,its,ids,uts)

      USE mo_bas_gsm, ONLY: pi,om,bk,re,amO2
      dimension ano2(its,ids,nh),tem(its,ids,nh)
     *         ,g(nh),rads(nh)
     *         ,qdis(2,its,ids,nh),gkoor(2,its,ids)
 
      dimension solu(nsu),sp(12)
      data 
cc      cross sections F10.7= 70
cc   *     sp/0.12e-18,1.39e-18,13.2e-18,10.5e-18,2.39e-18,
cc   *        0.87e-18,0.28e-18,0.01e-18,0.,0.,0.,0./
c       cross sections F10.7=115
     *     sp/7.29E-19,4.30E-18,4.03E-19,4.69E-19,2.29E-18,9.40E-18,
     *        1.37E-17,1.01E-17,5.98E-18,2.55E-18,1.08E-18,3.93E-19/
c       cross sections F10.7=180
cc   *     sp/8.14e-19,4.27e-18,4.51e-19,4.69e-19,2.31e-18,9.48e-18,
cc   *        1.37e-17,1.01e-17,6.01e-18,2.58e-18,1.09e-18,3.92e-19/
! cross section and spectral interval from Ackerman et al. Planet Space Sci. 1970. v. 1970 (20 intervals)
! and reduction to 12 intervals. Cross-sections averaging with flux value
!
! 1210.0  1220.0  
! 1220.0  1250.0
! 1250.0  1270.0
! 1270.0  1310.0
! 1310.0  1350.0  
! 1350.0  1380.0
! 1380.0  1500.0 
! 1500.0  1550.0
! 1550.0  1630.0 
! 1630.0  1670.0
! 1670.0  1720.0
! 1720.0  1760.0

      cr=pi/180.
      nsu05=nsu/2
      sum=0.
      do i=2,its-1
         do j=1,ids

	     gshir=gkoor(1,i,j)*cr
           gdol=gkoor(2,i,j)*cr
           gshir=0.5*pi-gshir
c     !!!!!!!  zenith angle   !!!!!
           coshi=sin(gshir)*sin(del)+cos(gshir)*cos(del)*
     *           cos(om*(uts-43200.)+gdol)
           hi=acos(coshi)
           sumI=0.
           do k=nh,1,-1
              ra=sqrt(rads(k)*(rads(k)+2.*re))
              alfa=atan(re/ra)
              him=pi-alfa
              if(hi.le.him) then
      
                hO2=(bk*tem(i,j,k))/(amO2*g(k))
                reh=(rads(k)+re)/hO2
              !  Chepmen function
                chep=chept(reh,hi)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                IF(K.EQ.NH) THEN 
                  sumI=anO2(i,j,k)*hO2
                else
                  sumI=sumI+0.5*(anO2(i,j,k)+anO2(i,j,k+1))*
     *                 (rads(k+1)-rads(k))    
                ! print sumI 
                end if
                sumL=0. 
                sumErg=0.                                               
                do l=1,nsu05                                   
                   tau=sp(l)*sumI*chep 
                   expTAU=exp(-tau)                                  
                   sumL=sumL+sp(l)*solu(l)*expTAU 
                   sumErg=sumErg+sp(l)*solu(l+nsu05)*expTAU                
                    ! print*,tau,k,l  
                end do                                  
                qdis(1,i,j,k)=sumL*1.e9 ! *anO2(i,j,k)!  
                qdis(2,i,j,k)=sumErg *anO2(i,j,k) 
              end if
            end do                           
          end do
      end do
      
      return 
      end

