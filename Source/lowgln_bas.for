c . . . ver. 2012 - P_riliv 
      subroutine lowgln_bas(pgl,rads,kpars,nh,its,ids,day
     *           ,ap,fa,fs,gkoor,dtets,ddolgs,uts,musl,pril,KPA,NT)
c
c     . . . parameters at 80 km altitude
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
        pgl(3,1,i,j)=5.4e+9
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