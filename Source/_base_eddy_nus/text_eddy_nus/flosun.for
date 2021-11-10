!
! model Nusinov based on paper 3.	Nusinov A. A. Models for prediction of 
! EUV- and X-ray solar radiation based on 10.7 cm radio emission // Proceedings 
! of the Workshop on the Solar Electromagnetic Radiation Study for Solar Cycle 
! 22 edited by R.F. Donnely
!
      subroutine flosuN(ff,fa0,sole,mas4)
      dimension al(35),bl(35),sole(mas4),bla(3),blr(3),
     *          sol(38),soler(8),ale(8),d(8),solei(8)

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

!  
!   ff -  bakground F10.7
!   fa0 - daily F10.7  
! 
!   F10.7 background  min = 63
      if(ff.lt.63) ff=63.
c
      fff=(ff-60.)**2
      fff=fff**0.333
      ffa=(fa0-ff)**2
      ffa=ffa**0.333
      fir=blr(1)+blr(2)*fff+blr(3)*ffa
!
!  base Nusinov 35th spectral interval
!   lam1 (A)     lam2 (A)
! 1   100       150    
! 2   150       200    
! 3   200       250    
! 4   256.32    256.32 
! 5   284.15    284.15 
! 6   250       300    
! 7   303.78    303.78 
! 8   300       350    
! 9   368.07    368.07 
! 10  350       400    
! 11  400       450    
! 12  465.22    465.22 
! 13  450       500    
! 14  500       550    
! 15  554.37    554.37 
! 16  584.33    584.33 
! 17  550       600    
! 18  609.76    609.76 
! 19  629.73    629.73 
! 20  600       650    
! 21  650       700    
! 22  703.36    703.36 
! 23  700       750    
! 24  765.15    765.15 
! 25  770.41    770.41 
! 26  789.36    789.36 
! 27  750       800    
! 28  800       850    
! 29  850       900    
! 30  900       950    
! 31  977.02    977.02 
! 32  950       1000   
! 33  1025.72   1025.72
! 34  1031.91   1031.91
! 35  1000      1050   
!
      do i=1,35
        sol(i+2)=(al(i)+bl(i)*fir)*fir
      end do 

!  Layman Alpha model based on Nusinov A.A., Katyushina V.V. 
!  Lyman alpha line intensity as a solar activity index. Solar Phys.
!  1994. V. 152. V 1. P. 201-206. 
      sol(38)=bla(1)+bla(2)*fff+bla(3)*ffa
!
!  XUV-model based on Nusinov
!
      ri820=0.023*fa0-1.440
!
      do i=1,8
        d(i)=1.56/ale(i)+0.22
        solei(i)=soler(i)*(ri820/soler(1))**d(i)*1.e-2
      end do
! 2 first XUV intervals
      sol(1)=solei(1)+solei(2)+solei(3)
      sol(2)=0.
      do i=4,8
        sol(2)=sol(2)+solei(i)
      end do
c
!  Combine 35 band sol to 15 band sole 
      sole(1)=sol(1)
      sole(2)=sol(2)+sol(3)/2.
      sole(3)=sol(3)/2.+sol(4)/3.
      sole(4)=sol(4)/3.
      sole(5)=sol(4)/2.
      sole(6)=sol(5)+sol(6)+sol(8)/2.
      sole(7)=sol(7)+sol(8)/2.+sol(9)+sol(10)/4.
      sole(8)=sol(10)*3./4.+sol(11)+sol(12)+sol(13)+sol(15)/4.
      sole(9)=sol(14)+sol(15)*3./4.
      do i=16,21
        sole(9)=sole(9)+sol(i)
      end do
      sole(9)=sole(9)+sol(22)/2.
      sole(10)=sol(22)/2.
      do i=23,29
        sole(10)=sole(10)+sol(i)
      end do 
      sole(11)=sol(30)+sol(31)+sol(32)/4.
      sole(12)=sol(32)*3./4.
c     do 6 i=33,37
c       sole(12)=sole(12)+sol(i)
c   6 continue
      sole(12)=sol(36)+sol(37)
      sole(13)=sol(33)
      sole(14)=sol(35)
      sole(15)=sol(38)
      return
      end
