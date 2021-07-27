      subroutine flosuEUVAC(ff,fa0,sole,mas4)
! 37 wave interval from EUVAC to 15 GSM
      dimension al(35),bl(35),sole(mas4),bla(3),blr(3),
     *          solEUVAC(37),sol(37),soler(8),ale(8),d(8),solei(8)
      integer godi
      data pi/3.1415926/,gb/1965.8/,t/10.2/
      data al /0.1144,1.6558,1.186 ,0.3672,-0.3421, 0.2435,
     *         5.8248,0.1965,0.5578,0.1043, 0.2446, 0.243,
     *         0.3089,0.3837,0.5837,1.000,  0.294,  0.4431,
     *         1.266, 0.2320,0.1571,0.313,  0.1125, 0.1478,
     *         0.133, 0.6098,0.5095,1.312,  3.037,  2.377,
     *         3.847, 1.071, 2.944, 1.633,  1.5767/
      data bl/-0.0018,-0.0932, 0.0879,-0.0289, 0.389,   0.5938,
     *        -0.3854, 0.321, -0.0487, 0.1061, 0.02153,-0.0176,
     *         0.0658,-0.0039,-0.0592, 0.,    -0.012,   0.0283,
     *        -0.0876, 0.0183,-0.0067,-0.0342,-0.0003, -0.0101,
     *         0.0288,-0.0631, 0.0109,-0.046, -0.0404, -0.022,
     *        -0.353,-0.03296,-0.0305,-0.065, -0.0471/
      data bla/134.,39.9,5.47/,blr/0.725,0.160,0.0592/
      data soler/1.32,7.4,5.6,12.5,12.9,15.6,17.3,18.1/,
     *       ale/20.,40.,50.,60.,70.,80.,90.,100./

c . . . F10.7 - минимальный = 63 
      if(ff.lt.63) ff=63.  ! background f107
c
      fff=(ff-60.)**2
      fff=fff**0.333
      ffa=(fa0-ff)**2      ! fa0 - current f107
      ffa=ffa**0.333
      fir=blr(1)+blr(2)*fff+blr(3)*ffa
c
!      do 1 i=1,35
!        sol(i+2)=(al(i)+bl(i)*fir)*fir  !!! 37 bins spectra
!    1 continue
      call euvac(ffa0,ff,solEUVAC)

      solLA=bla(1)+bla(2)*fff+bla(3)*ffa !!! Layman Alpha ???
c
      ri820=0.023*fa0-1.440
      do 2 i=1,8
        d(i)=1.56/ale(i)+0.22
        solei(i)=soler(i)*(ri820/soler(1))**d(i)*1.e-2
    2 continue
      sol(1)=solei(1)+solei(2)+solei(3)
!      sol(2)=0.
      do 3 i=4,8
        sol(2)=solEUVAC(1)+solei(i)
    3 continue
!     correction wave interval according to FLOSU
      do i=3,8
        sol(i)=solEUVAC(i-1)*2.0 ! increase for Nusinov value
      end do
      do i=9,37     
        sol(i)=solEUVAC(i)*2.0   ! increase for Nusinov value 
      end do
c
      sole(1)=sol(1)
      sole(2)=sol(2)+sol(3)/2.
      sole(3)=sol(3)/2.+sol(4)/3.
      sole(4)=sol(4)/3.
      sole(5)=sol(4)/2.
      sole(6)=sol(5)+sol(6)+sol(8)/2.
      sole(7)=sol(7)+sol(8)/2.+sol(9)+sol(10)/4.
      sole(8)=sol(10)*3./4.+sol(11)+sol(12)+sol(13)+sol(15)/4.
      sole(9)=sol(14)+sol(15)*3./4.
      do 4 i=16,21
        sole(9)=sole(9)+sol(i)
    4 continue
      sole(9)=sole(9)+sol(22)/2.
      sole(10)=sol(22)/2.
      do 5 i=23,29
        sole(10)=sole(10)+sol(i)
    5 continue
      sole(11)=sol(30)+sol(31)+sol(32)/4.
      sole(12)=sol(32)*3./4.
c     do 6 i=33,37
c       sole(12)=sole(12)+sol(i)
c   6 continue
      sole(12)=sol(36)+sol(37)
      sole(13)=sol(33)
      sole(14)=sol(35)
      sole(15)=solLA
      return
      end

      SUBROUTINE EUVAC(F107,F107A,EUVFLX)
      INTEGER I
      REAL F107,F107A,EUVFLX(37),AFAC(37),F74113(37),FLXFAC
C
C------ F74113 reference spectrum (doubled below 150-250 A, tripled <150)
C------ Will be multiplied by 1.0E9 later
      DATA F74113/1.20,0.450,4.800,3.100,0.460,0.210,1.679,0.8
     > ,6.900,0.965,0.650,0.314,0.383,0.290,0.285,0.452,0.720
     > ,1.270,0.357,0.530,1.590,0.342,0.230,0.360,0.141,0.170
     > ,0.260,0.702,0.758,1.625,3.537,3.000,4.400,1.475,3.500
     > ,2.100,2.467/
C
C--- Scaling factors(Ai) for the EUV flux
      DATA AFAC/1.0017E-02,7.1250E-03,1.3375E-02,1.9450E-02,2.7750E-03
     > ,1.3768E-01,2.6467E-02,2.5000E-02,3.3333E-03,2.2450E-02
     > ,6.5917E-03,3.6542E-02,7.4083E-03,7.4917E-03,2.0225E-02
     > ,8.7583E-03,3.2667E-03,5.1583E-03,3.6583E-03,1.6175E-02
     > ,3.3250E-03,1.1800E-02,4.2667E-03,3.0417E-03,4.7500E-03
     > ,3.8500E-03,1.2808E-02,3.2750E-03,4.7667E-03,4.8167E-03
     > ,5.6750E-03,4.9833E-03,3.9417E-03,4.4167E-03,5.1833E-03
     > ,5.2833E-03,4.3750E-03/
C
C----- loop through the wavelengths calculating the scaling factors and
C----- the resulting solar flux.
C----- The scaling factors are restricted to be greater than 0.8
       DO 50 I=1,37
          FLXFAC=(1.0 + AFAC(I) * (0.5*(F107+F107A) - 80.0))
          IF(FLXFAC.LT.0.8) FLXFAC=0.8
 !         EUVFLX(I)=F74113(I) * FLXFAC * 1.0E9
          EUVFLX(I)=F74113(I) * FLXFAC
50    CONTINUE
      RETURN
      END