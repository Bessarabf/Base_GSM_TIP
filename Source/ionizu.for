      subroutine ionizu(sole,solen,rads,par,gkoor,ut,dolm,
     *       ddolgs,del,nh,its,ids,nse,kpars)
!      implicit none
!*** Start of declarations inserted by SPAG
!      real A1 , A2 , A3 , AB , AC , AE , AI , AM , BK , BKG , CHAB , 
!     &     CHEP , CHEPT , CONE , CR , CX , DDOLGS , DEL , DOLGG , DOLM
!      real EX1 , EX2 , FS , G0 , GKOOR , HL , HSM , OM , PAR , PI , 
!     &     RADS , RE , REH , RLAT , SA , SI , SM , SOLE , SOLEN , SX
!      real UT , XI , XIG , XLG
!      integer I , IDS , IG , II , ITS , IV , J , K , KL , KPARS , L , 
!     &        M , MM , ND , NH , NSE
C*** End of declarations inserted by SPAG
      dimension am(6) , rads(nh) , ai(4) , sa(4,15) , si(4,15) , 
     &          sole(nse) , solen(nse) , gkoor(2,its,ids) , 
     &          par(kpars,nh,its)
c     double precision tau0,xp,xxp,sx,reh,chep,derf,hsm,
c    *       re,pi,xlg,cr,xig,xi,cx,chab,ab,ex1,ex2
      REAL*8 ex1, ex2, fs ! ab, chab, 
      REAL*8,PARAMETER :: bk = 1.38D-16, ae=1.66D-24, g0 = 980.665D0,
     &  pi = 3.1415926D0, re = 6371.02D5,  om = 7.2722D-5
      data am/32. , 28. , 16. , 30. , 4. , 1./ 
!      , bk/1.38E-16/ , 
      data si/4.63, 4.98, 2.49, 0.0,  7.0, 4.23, 6.0, 0.0 ,12.68, 9.4,
     *        7.20, 0.0 ,13.0 , 9.0,  7.5, 0.0 ,12.5, 9.36, 7.58, 0.0,
     *        13.0, 8.56, 9.12, 0.0, 15.5, 8.5 , 9.2, 0.0 ,19.53,18.00,
     *        10.5, 0.0 ,24.12,23.47,12.69,0.0 ,16.27,20.58,
     *        8.16, 0.0 , 5.93, 0.0 , 4.0 ,0.0 , 4.8,3*0.0 ,2.5,0.0,
     *        2*0.0,1.,0.0,5*0.,2.02/
      data sa/0.46, 0.34, 0.16, 0. , 1.67, 0.8, 1.00, 0. , 4.68 ,3.48,
     *        4.24, 0.  , 6.5 , 5.0, 5.0 , 0. , 8.5 , 5.4, 6.50, 0.0 ,
     *        12.0, 7.20, 8.00, 0.0, 15.5, 8.5, 9.2 , 0. ,19.00,18.50,
     *        10.5, 0.  ,24. , 24.4, 12.7, 0. ,27.3 ,27.3, 8.50, 0.0 ,
     *        16.2, 9.0 ,  4.0, 0. , 7.7 , 0.02, 0. , 0.0, 4.00, 0.0 ,
     *         0. , 0.  ,  1.6, 0.0, 0.0 , 0.0 , 0.01,0.0, 0.0 , 2.42/

      cr = 180./pi
      bkg = bk/g0/ae
      k = 4
!      kl = nse
      kl = 15
      nd = dolm/ddolgs + 1

      do ig = 1 , its
!        PRINT *,'par(1) max=',maxval(par(1,:,ig)),'min=',minval(par(1,:,ig))
!        PRINT *,'par(2) max=',maxval(par(2,:,ig)),'min=',minval(par(2,:,ig))
!        PRINT *,'par(3) max=',maxval(par(3,:,ig)),'min=',minval(par(3,:,ig))
!        PRINT *,'par(4) max=',maxval(par(4,:,ig)),'min=',minval(par(4,:,ig))
!        rlat = gkoor(1,ig,nd)/cr
!        rlat = pi/2. - rlat
        rlat = pi/2. - gkoor(1,ig,nd)/cr
!        dolgg = gkoor(2,ig,nd)
!        hl = ut + dolgg*3600./15.
        hl = ut + gkoor(2,ig,nd) * 240 ! 3600. / 15. = 240
        cx = sin(del)*sin(rlat) + cos(del)*cos(rlat)*cos(om*(hl-43200.))
!        sx = sqrt(1.0-cx*cx)
        xi = acos(cx)
c       if(cx.lt.0.) xi=pi-xi
        xig = xi*cr
c ***** for Layman night radiation only
        if ( xig>90.0 ) then
          xlg = -0.96*cos(1.2*(xig-180.0)/cr)
        else
          xlg = 0.0
        endif
!        PRINT *,'hl=',hl,'om*(hl-43200.)=',om*(hl-43200.),'xlg=',xlg
!------------------------------------------------------------------------------
!        do i = 1 , k
          ai(:) = 0.
!        enddo
        do iv = 1 , nh
          i = nh - iv + 1
          cone = 0.;  sm = 0.
          do j = 1 , k
            cone = cone + par(j,i,ig)
            sm = sm + am(j)*par(j,i,ig)
          enddo
          sm = sm/cone
!          hsm = sm/(bkg*par(7,i,ig))
!          reh = (rads(i)+re)*hsm
          reh = (rads(i)+re)*sm/(bkg*par(7,i,ig))
c ********* Chepmen function**********
          if ( xig<80.0 ) then
            chep = 1.0/cx
          ELSE
            chep = chept(reh,xi)
          endif
!          if ( i/=nh ) then
          if ( i/=nh ) ai(1:k) = ai(1:k) + (par(1:k,i,ig) + 
     &         par(1:k,i+1,ig)) * (rads(i+1)-rads(i)) * 0.5
!	      if(ai(m).lt.0) print*, ai(m),m,par(m,i,ig),i,ig
          do j = 1 , k
            fs = 0.D0
            do l = 1 , kl

              ab = 0.
              do m = 1 , k
                ab = ab + sa(m,l)*ai(m)
              enddo
              ab = ab*1.E-18
              chab = chep*ab
             if(chab >= 10.0) chab = 10.0  !chab=100.0
              ex1 = exp(-chab)
             if(ab >= (10.0 + xlg) ) then  !ab=100.0
                if (xlg>=0) then
                  ab = 10.0 
                else 
                  ab = 10.0 + xlg
                endif
              endif
              ex2 = exp(-ab + xlg )
!              if ( ex1<=1.E-30 ) ex1 = 1.E-30
!              if ( ex2<=1.E-30 ) ex2 = 1.E-30
              fs = fs + DBLE(si(j,l))*(DBLE(sole(l))*ex1 + DBLE(solen(l))*ex2)
            enddo ! l
            fs = fs*1.D-9
            mm = j + 12
            if ( j==3 ) mm = 16
            if ( j==4 ) mm = 15
            par(mm,i,ig) = par(j,i,ig) * REAL( fs )
            if ( par(mm,i,ig)<=1.E-20 )  par(mm,i,ig) = 1.E-20 ! 1.E-6
!            if (par(mm,i,ig) > 20000.) then
!              PRINT *, i,ig
!              stop
!            endif
          enddo ! j
        enddo ! iv
!        if ( chab > 10)  PRINT *,'ex1=',ex1,'ex2=',ex2,
!     & 'fs=',fs,'chab=',chab,'ab=',ab,'chep=',chep,'ai=',ai 
        do ii = 13 , 16
          a1 = alog(par(ii,nh-1,ig))
          a2 = alog(par(ii,nh-2,ig))
          a3 = alog(par(ii,nh-3,ig))
          ac = DBLE(a3) + 3.D0 * DBLE(a1) - 3.D0 * DBLE(a2)
          par(ii,nh,ig) = exp( REAL(ac) )
!          if (par(ii,nh,ig) > 20000.) PRINT *,'a1=',a1,'a2=',a2,'a3=',a3,'ac=',ac
!          if (par(ii,nh,ig) > 20000.) PRINT *,'a1=',par(ii,nh-1,ig),
!     & 'a2=',par(ii,nh-2,ig),'a3=',par(ii,nh-3,ig),'ac=',exp(ac)
!          if (ii == 14) PRINT *,'a1=',a1,'a2=',a2,'a3=',a3,'ac=',ac
!          if (ii == 14) PRINT *,'a1=',par(ii,nh-1,ig),
!     & 'a2=',par(ii,nh-2,ig),'a3=',par(ii,nh-3,ig),'ac=',exp(ac)
        enddo ! ii

      enddo ! ig
      end subroutine IONIZU