c   terpot_HAM - integrate with HAMMONIA 
c   ver.    18.05.18 Ion Drag sent to HAMMONIA 
c   ver.    10.05.18 Joul heating sent to HAMMONIA
c   ver.    08.04.14 allocatable massives
c   version 25.05.12 add to intrface KPA & NT for massive pril
c           20.11.17 integate with HAMMONIA
c - terpot
c - lowgln
c - gstrf0
c - pgl3d
c - timol
c - nachc
c - plots
c - turbk
c - conn
c - cntn90
c - conno1
c - co2con
c - rezam
c - tnalt
c - nts
c - bongl
c - boskli
c - progjn
c - sumro
      subroutine terpot_ham(day,god,dayt,godt,uts,tau,dts,solet,sole,
     *       solu,nsu,nse,kpars,rads,nh,gkoor,its,ddolgs,dtets,fa,fs,
     *       ap,pkp,dst,ae,al,au,bmpz,bmpy,mass,delta,pgl,pgi,ids,ins
     *      ,isp,vir,vid,vim,verno,parj,potef,ntr,nl2,pril,KPA,NT)
!    
!     module with mass GSM to HAM & HAM to GSM 
      USE mo_ham_gsm

!	  use ieee_arithmetic
!	  Parameter(nh0=30,its0=37,ids0=72) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
     *   ,q(:,:,:),g(:),rp(:),ctd(:),anco2(:,:,:)
     *   ,cHot(:,:,:),tHot(:,:,:),ron(:,:,:),ros(:,:,:),ros0(:,:,:)
    
       allocate (an1(its,ids,nh),an2(its,ids,nh),an3(its,ids,nh) 
     *   ,an6(its,ids,nh)
     *   ,an11(its,ids,nh),an21(its,ids,nh),an31(its,ids,nh)
     *   ,an61(its,ids,nh),vr(its,ids,nh),vi(its,ids,nh)
     *   ,vj(its,ids,nh),vi1(its,ids,nh),vj1(its,ids,nh)
     *   ,q(its,ids,nh),g(NH),rp(NH),ctd(NH),anco2(its,ids,nh)
     *   ,cHot(its,ids,nh),tHot(its,ids,nh)
     *   ,ros0(ITS,IDS,NH),ros(ITS,IDS,NH),ron(ITS,IDS,NH))

      data pi,om/3.1415926,7.27e-5/,! bk/1.38e-16/,
     *     g0,re/981.,6.371e8/
     *    ,key/1/

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
!      rate dissociation massiv
ccc   	call r_dis(qdis,pgl,solu,
ccc     *           gkoor,rads,delta,kpars,nh,its,ids,nsu,uts)
       
       call lowgln_ham(pgl,rads,kpars,nh,its,ids,day
     *         ,ap,fa,fs,gkoor,dtets,ddolgs,uts,mass(18),pril,KPA,NT)

       call pgl3d(pgl,kpars,nh,its,ids,an1,an2,an3,an6,vr,vi,vj)
   
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
!!!      Fist step on time
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
        call nonew(pgl,pgi,cHot,tHot,gkoor,ctd,rads,rp,g,
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
         call heatpo_ham(pgl,pgi,parj,solet,solu,nsu,nse,
     *               kpars,rads,g,nh,gkoor,its,ddolgs,dtets,
     *               mass,delta,day,uts,tau,dts,ctd,vim,vid,
     *               vir,ids,ins,an1,an2,an3,anco2,an6,an61,ros,
     *               vi,vj,vr)


!  correction Tn by HAMMONIA
      ncor=10   !   10 ! point number of correction 
      do i=1,its
	  do j=1,ids
	    do k=1,ncor-1
		  an61(i,j,k)=gsmHAM(k,i,j)
	    end do
            an61(i,j,ncor)=(gsmHAM(ncor,i,j)+
     *                         an61(i,j,ncor))*0.5
	   end do
	end do
!  correction end
         call nts (an61,nh,its,ids,nh,its-3)

      end if

      if(mass(5).ne.0) then
         call sdizkn_ham(an1,an2,an3,an11,an21,an31,
     *               an61,vr,vi,vj,ros,rp,rads,g,nh,its,ids,
     *               dts,ctd,roS,solu,gkoor,delta,nsu,dtets,uts,ddolgs)
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
! correction Vn 
         do i=1,its
           do j=1,ids
             do k=1,ncor-1
               vi1(i,j,k)=UgsmHAM(k,i,j)
               vj1(i,j,k)=VgsmHAM(k,i,j)
             end do
             vi1(i,j,ncor)=(UgsmHAM(ncor,i,j)+
     *                         vi1(i,j,ncor))*0.5
             vj1(i,j,ncor)=(VgsmHAM(ncor,i,j)+
     *                         vj1(i,j,ncor))*0.5
           end do
         end do
!  correction end
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
     *   ,q,g,rp,ctd,anco2
     *   ,cHot,tHot,ron,ros,ros0)

      return
      end
c
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
c
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
c
      subroutine conn(pgl,kpars,nh,its,ids)
       dimension pgl(kpars,nh,its,ids)
          do 1 k = 1 , nh
          do 1 i = 1 , its
          do 1 j = 1 , ids
       pgl(5,k,i,j)=0.
  1     continue
       return
       end



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
c

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
c
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
c
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
c
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
          s=alog(ro(i,j,k-1)*(tn/tv)*(amcv/amcn))-alf
          ro(i,j,k)=exp(s)
    3    continue
    2  continue
    1 continue
      return
      end
c
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

	subroutine bospgl(pgl,kpars,nh,its,ids,np)
c . . . Расчет в полюсной точке. Посвящается Клименко В.И.
C   . . . ids - четное
      dimension pgl(kpars,nh,its,ids)
      i2=its-1
      do 1 k=1,nh
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
c
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
c
      subroutine turbk(ctd,rads,n)
      dimension ctd(n),rads(n)
c      data hm/103./,c1,c0/1.e 6,5.e 5/,
c      data hm/ 94./,c1,c0/1.e 6,5.e 5/,
c    *     s1,s2,s3/0.05,0.05,0.07/
c      . . . Вариант с увеличенным Кт и медленным спадом
c  !     data hm/ 94./,c1,c0/5.e 6,5.e 6/,
c      data hm/ 94./,c1,c0/1.e 6,5.e 5/,
      data hm/ 90./,c1,c0/5.e 6,1.e 6/,!9.03.11
     *     s1,s2,s3/0.01,0.05,0.01/
!	data hm/98./,c1,c0/5.0e 6,3.0e 6/
!     *     s1,s2,s3/0.01,0.005,0.01/
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
c
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
c
      subroutine cntn90(rads,fig,dolg,td,ap,fa,fs,tlttau,utsec,
     *                 cn1,cn2,cn4,tn,nh)
c   cn4,cn2,cn1,dh,tn - мaccивы концентраций O,O2,N2,H и тeмпepaтyp
c    в yзлax на высотах rads в точке fig,dolg в момент td,uttau
      dimension rads(nh),cn1(nh),cn2(nh),cn4(nh),tn(nh)
      dimension d(8),t(2),apm(7)
      
c     do 3 i=1,7
c       apm(i)=ap
c 3   continue
      apm(1)=11.
      apm(2)=2.
      apm(3)=4.
      apm(4)=6.
      apm(5)=12.
      apm(6)=13.
      apm(7)=44.
      alat=fig*1.7453292e-2
      alon=dolg*1.7453292e-2
      iyd=80*1000+td
      do 2 ih=1,nh
        call gtd6(iyd,utsec,rads(ih)/1.e5,fig,dolg,tlttau,fs,fa,apm
     *      ,48,d,t)
         tn(ih)=t(2)
         cn1(ih)=d(3)
         cn2(ih)=d(4)
         cn4(ih)=d(2)
 2    continue
      return
      end

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
