c                       CYCLT
C    22/05/2014 ADD PAR2 TO INTERFACE
c    16/04/14 add nl2 to interface
C    25/05/2012 - added kpa & nt - parameters of priL and AE
      SUBROUTINE CYCLT_BAS(god,day,ut0,utk,dtt,dts,tau,b,c,sole,solu,bmpz,bmpy,
     &                     bmod,vsol,csol,solen,solet,fa0,fs,ap0,pkp0,dst0,ae0,
     &                     al0,au0,nh,ddolgs,dtets,ddolgt,ntsl,nl,nl2,q,u,dtett,
     &                     int,ins,rads,ntr,gkoor,kpart,kpars,nadrt,nadrs,ns,
     &                     idt,ids,its,nzapt,nzaps,mas,mast,mass,uprt,kdf,kdu,
     &                     ldor,isp,rmaxt,par1,PAR2,par,pari,pole,nr,verno,ni,
     &                     park,ks,nv,nsu,nse,tet1,tet2,tet3,pdpc,eps0,om,fac1,
     &                     fac2,fac3,pglo,imja,dut,hut,iput,keut,koob,kut,nara,
     &                     naou,sut,tut,izap,pril,kpa,nt,nzapjet,lzapjet,isat)
!------------------------------------------------------------------------------
      ! god � ���, �����, ������ ���� 
      ! day � ���� ���� (�� 1-�� ������)
      ! ut0 � �����, ��� ������ �������
      ! utk � �����, ���, ��������� �������
      ! dtt � ��� �� �������, ���, � �������������� ��������� �������
      ! dts - ��� �� �������, ���, � �������������� ��������� �����
      ! tau - ��� �� �������, ���. (��������������, ��� tau, dtt,dts ����� ����������. ������ ��� ���� �����).
      ! b - 
      ! c - 
      ! sole � ������ ����������� nse (15 � ������� ��������) � ������� ������������� ��������� (���) � ������� (*1.0e11)
      ! solu � ������ ����������� 2*nsu (2*12) � ������� ��������������� ��������� ������ � ������� (nsu) � ����� (nsu).
      ! bmpz � bz-���������� �������������� ���������� ���� (���)
      ! bmpy � by-���������� ���
      ! bmod � ������ ���,
      ! vsol � �������� ���������� �����,
      ! csol - ,
      ! solen - ������ ����������� nse � ������ ����������� ���
      ! solet - ������ ����������� nse � �� ������������,
      ! fa0 � ������ F10.7 ������� (��������������),
      ! fs � ������ F10.7 ������� (��� ������ ��������) ��� �� 81 ���� ��� EUVAC. 
      ! ap0 � Ap ������,
      ! pkp0 � Kp ������,
      ! dst0 � Dst ������,
      ! ae0 � AE -������,
      ! al0 � Al ������,
      ! au0 � AU ������, 
      ! nh � ����� ����� �� ��������� � ����� (30),
      ! ddolgs � ��� �� ������� � �����, ������
      ! dtets � ��� �� ������ � �����, ������
      ! ddolgt � ��� �� ������� � �������, ����,
      ! ntsl - ����� ����� � ������� ������
      ! nl - ����� ������� ����� � ��������� ����� ��������� ��������� (17 ��� 5 ���� ����� �� ������),
      ! nl2 - nl+nl+3 ,
      ! q - ������ q(nv,nl)-��������� ��������� ������� ��������� (����� ������� �����, ���� ����� ������� �����)
      ! u - ������ u(nl) - ��������� ��������� ������� ��������� ���������������� ������� �����,
      ! dtett -  � ��� �� ������ � �������, ����,
      ! int - ����� ���������������� ���������� � ������� �� ������(����� ,
      ! ins - ����� ���������������� ���������� � ����� �� ������� (����� 6),
      ! rads - ������ �������� ����� ����� ������������� ����� ������, ����� ��������� nh=30 , �� 80 ��  526 �� (� �����������)
      ! ntr - ����� ����� ����� ������� ����� ��������� ������� ���������,
      ! gkoor - ������ ����������� (2, its,ids) ������������ ����� ������������ ��������� (�������� � �������, ����) � ��������������.
      ! kpart - ����� �������������� ���������� � ������� ,
      ! kpars - ����� �������������� ���������� � ����� (19),
      ! nadrt - ����� ��������� ������ ���������� ������� � f4 ,
      ! nadrs - ����� ��������� ������ ���������� ������ � f4 ,
      ! ns - ,
      ! idt - ����� ��������� ����� � ������� (24),
      ! ids - ����� ��������� ����� � ����� (24),
      ! its - ����� �������� ����� � ���� (37) ,
      ! nzapt - ,
      ! nzaps - ,
      ! mas - ������ ����������� ���������� ������ (10),
      ! mast - ������ ����������� ���������� ��� ��������� ������� ������� (40),
      ! mass - ������ ����������� ���������� ��� ��������� ������� ������ (30),
      ! uprt,
      ! kdf - ������ ������� ������� � f4, ��������������� ������ ���������� ������, ����������, ������� � ��.  
      ! kdu - ������ ����� ������� � f4,  ��������������� ���������� ������, ����������, ������� � ��.  
      ! ldor - ����� ������ � f4 � ������ (4096),
      ! isp,
      ! rmaxt,
      ! par1 - ������ ���������� ���� � ����� ��������� ��������� (��������, ��������),
      ! PAR2 - ������ ���������� ���� � ����� ��������� ���������,
      ! par - ������ ���������� ���� � ����� ��������� ���������,
      ! pari - ������ ����������� ni - ��������� ������ � ����� ��������� ���������,
      ! pole - ������ ������������ 1024 ��� ������ ���������� �� ����� f4,
      ! nr - ������ ������� ���������� ���� (���� ��������� ���������) =  nh*its*kpars ,
      ! verno - ��������������� ���������� ���������� (.TRUE. ��� .FALSE.),
      ! ni - ������ ������� ������ ( ���� ��������� ���������)  ������ ntpl*int ,
      ! park,
      ! ks - ����� ����� � ��������� ��������� ������� ntpl*2,
      ! nv,
      ! nsu - ����� ������������ ���������� ��� (12),
      ! nse - ����� ������������ ���������� ��� (15 � ������� ��������),
      ! tet1 - �������� ��������� ���������� ����� 1,
      ! tet2 - �������� ��������� ���������� ����� ���� 2,
      ! tet3 - �������� ��������� ���������� ����� ���� 3 (���� ����)
      ! pdpc - �������� ����������� ����� ����� � ��,
      ! eps0,
      ! om,
      ! fac1 - �������� ����� ���� 1,
      ! fac2 - �������� ����� ���� 2,
      ! fac3 - �������� ����� ���� 3,
      ! pglo - 4-� ������ ������, ���������� ������, ����������� (kpars,nh,its,ids)
      ! imja,
      ! dut - ������� �����, ��������� ��� ������ ���������� (��������������� ����������),
      ! hut - ������ �����  ��������� ��� ������ ���������� (��������������� ����������),
      ! iput ,
      ! keut - ����� ������ � ���� f4 ,
      ! koob - ����� ������� (?), ��������� ��� ������,
      ! kut - 
      ! nara -��� �������� ����� � ������������ ������� ,
      ! naou - ��� ��������� ����� � ������������ �������,
      ! sut - ������ ������� (?),
      ! tut - ����� �������,
      ! izap ,
      ! pril - ������ ���������� ��� ��������� ������� ������ (������� �� ������ �������) , ����������� kpa,its,ids,nt
      ! kpa - ����� ���������� ��������� ������� ,
      ! nt - ����� ��������� ������� ��� pril, ������ 12 ��� 24
      ! nzapjet,
      ! lzapjet,
      ! isat - 
!------------------------------------------------------------------------------
      CHARACTER imja*80
      CHARACTER nara(100)*25, naou(100)*25
      INTEGER god, day, uprt, verno, dayt, godt
      LOGICAL readfl
      DIMENSION sole(nse), solu(nsu), solen(nse), solet(nse), mas(10), mast(40), 
     &          mass(30), ps(10), u(nl), kdf(20), kdu(20), pole(ldor/4), par(nr)
     &          , par1(nr), par2(nr), ntsl(nl), q(nv,nl), rads(nh), 
     &          dut(isat,100), hut(isat,100), iput(100), keut(100), kut(100), 
     &          izap(100), sut(isat,100), tut(isat,100), nzapjet(5), lzapjet(5)
      DIMENSION gkoor(2,its,ids), pari(ni), park(ks), pril(*)
      DIMENSION pglo(kpars,nh,its,ids)
!
      ALLOCATABLE vdr(:), potef(:,:,:), gins(:,:,:,:)
      ALLOCATE(vdr(ks),potef(ntr,idt,nl2),gins(ins,nh,its,ids))
      PRINT *, 'cyclt - begin'
 
      readfl = .TRUE.
      nfile = 5
      CALL ZAIT(readfl,nfile,nzapt,nzaps,nadrt,nadrs,isp,ldor,kdf,god,day,ut0,
     &          ut1,pole)
      uut1 = ut1 + 0.0001
      IF( uut1>utk )THEN
          PRINT 171, ut1, utk
 171      FORMAT('GSMTIP: cyclt - ut1='g10.2,' > utk=',g10.2/
     &           '***** incorrect start and stop time of calculating! ***** ')
          STOP
      ENDIF ! uut1>utk
 !     nl2=nl*2+3
      IF( mast(23)==0 )THEN
! . . . without dreif
          potef = 0.     ! potef(kkk,jjj,iii)=0.
      ELSE
! . . . potential calculate
          readfl = .TRUE.
          n3 = ntr * idt * nl2
          CALL WPOTEF(readfl,potef,n3,kdf,kdu,ldor,isp)
      ENDIF ! mast(23)==0
      iqo = mast(25)
      ut = ut0
      nfr = 6
      nfw = 13
! . . . reading init cond and writing in iteratioon part of tube
      CALL COPMD(nfr,nfw,kdf,isp,ldor,kdu,par,nr,mass,mast)
      vdr(1:ks) = 0.
      nomsl = 1
      md = 1
      dolm = 0.
      readfl = .TRUE.
      kpar = kpart
      nfile = 6
! . . . reading init cond and writing interpolation parameters
      DO i = 1, idt
          CALL WWT(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,kdf,ldor,isp,md,
     &             par,pole,nr,mast)
          CALL INTS(dolm,par,nr,rads,nh,ni,pari,ins,its,park,ks,ntsl,nl,ntr,kdf,
     &              ldor,isp,pole,kpart,vdr,dtett,u,ddolgs,dtets,gins,ids)
          dolm = dolm + ddolgt
      ENDDO
! . . . data in poles
      CALL BONGI(gins,ins,nh,its,ids)
!     . . . continue or new calculation?
      IF( uprt==3 )THEN
          CALL DAT(ut,dayt,godt,verno,ut0,ut1,day,god)
      ELSE
          ut1 = ut0
          dayt = day
          godt = god
! . . . date and time write in inf record
          READ(4,REC=1)pole
          pole(2) = god
          pole(3) = day
          pole(4) = ut0
          pole(5) = ut1
          WRITE(4,REC=1)pole
! . . .
      ENDIF ! uprt==3 
      delta = DEL(dayt)
      utn = ut1 + tau
      md = 0
      dtt0 = dtt
      IF( dtt0>dts ) dtt0 = dts
!!!!!!!!!!!! cycle on time begin  !!!!!!!!!!!!!!!!!!!!!
 121  IF( utn>=(utk+tau) ) GOTO 21
      uts = ut1 + dts
      IF( uts>=86400. )THEN
          CALL DAT(ut,dayt,godt,verno,ut0,uts,day,god)
          delta = DEL(dayt)
      ENDIF ! uts>=86400.
 122  IF( uts>=utn ) GOTO 22
!!!!!!!!!!!! from danmodel
      pkpj1 = pkp0
      AEpdpc = ae0
      AEj2 = ae0
!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL IACFLO_AE(uts,sole,solu,bmpz,bmpy,mas,csol,vsol,fa0,pkp0,ap0,ae0,
     &               dst0,al0,au0,fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps,AEj2)
      CALL CYCLT1_BAS(god,day,ut0,ut1,dtt,dts,tau,solen,sole,solu,solet,bmpz,
     &                bmpy,bmod,vsol,csol,ntr,gkoor,ps,gins,fa,fs,ap,pkp,dst,ae,
     &                al,au,rads,ni,nv,vdr,nh,ddolgs,dtets,ddolgt,ntsl,nl,idt,
     &                ids,its,q,u,kpart,kpars,nadrt,nadrs,nzapt,nzaps,mast,mass,
     &                mas,dtt0,uts,ns,dtett,int,ins,b,c,kdf,kdu,ldor,isp,rmaxt,
     &                nsu,nse,par1,PAR2,par,park,ks,pari,pole,nr,verno,delta,
     &                dayt,godt,tet1,tet2,tet3,pdpc,eps0,om,fac1,fac2,fac3,nl2,
     &                potef,iqo,pglo,imja,dut,hut,iput,keut,koob,kut,nara,naou,
     &                sut,tut,izap,pril,kpa,nt,nzapjet,lzapjet,AEpdpc,AEj2,
     &                pkpj1,isat)
      IF( verno==1 )uts = utn
      uts = uts + dts
      IF( uts>=86400. )THEN      ! ut hour > 24 h, new day
          CALL DAT(ut,dayt,godt,verno,ut0,uts,day,god)
          delta = DEL(dayt)
      ENDIF  ! uts>=86400.
      IF( uts>=utn )THEN
          dts = utn - uts + dts
          uts = utn
          IF( uts>=86400. )THEN
              CALL DAT(ut,dayt,godt,verno,ut0,uts,day,god)
              delta = DEL(dayt)
          ENDIF ! uts>=86400.
      ENDIF ! uts>=utn
      GOTO 122
 22   CONTINUE
!!!!!!!!!!!! from danmodel
      pkpj1 = pkp0
      AEpdpc = ae0
      AEj2 = ae0
!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL IACFLO_AE(uts,sole,solu,bmpz,bmpy,mas,csol,vsol,fa0,pkp0,ap0,ae0,
     &               dst0,al0,au0,fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps,AEj2)
      CALL CYCLT1_BAS(god,day,ut0,ut1,dtt,dts,tau,solen,sole,solu,solet,bmpz,
     &                bmpy,bmod,vsol,csol,ntr,gkoor,ps,gins,fa,fs,ap,pkp,dst,ae,
     &                al,au,rads,ni,nv,vdr,nh,ddolgs,dtets,ddolgt,ntsl,nl,idt,
     &                ids,its,q,u,kpart,kpars,nadrt,nadrs,nzapt,nzaps,mast,mass,
     &                mas,dtt0,uts,ns,dtett,int,ins,b,c,kdf,kdu,ldor,isp,rmaxt,
     &                nsu,nse,par1,PAR2,par,park,ks,pari,pole,nr,verno,delta,
     &                dayt,godt,tet1,tet2,tet3,pdpc,eps0,om,fac1,fac2,fac3,nl2,
     &                potef,iqo,pglo,imja,dut,hut,iput,keut,koob,kut,nara,naou,
     &                sut,tut,izap,pril,kpa,nt,nzapjet,lzapjet,AEpdpc,AEj2,
     &                pkpj1,isat)
      ut1 = utn
      utn = utn + tau
      readfl = .FALSE.
      nfile = 7
      nrazm = ins * nh * its * ids
      nkd = 20
      nf = 4
! . . . end time-step and writing into f4
      CALL READWR(readfl,nfile,kdf,kdu,nkd,nrazm,gins,ldor,isp,nf)
! . . . copy f4 or not?
      IF( mas(9)/=0 )THEN
          iuts = uts + 0.0001
          idts = mas(9)*3600 + 0.0001
          IF( MOD(iuts,idts)==0 )THEN
              isf4 = 0
              DO i = 1, 8
                  isf4 = isf4 + kdu(i)
              ENDDO
!            call copyf4(uts)
!            call copyf4day(uts,day)
              CALL CPF4DAY(uts,day,isf4,ldor)
          ENDIF ! MOD(iuts,idts)==0
      ENDIF !  mas(9)/=0
      GOTO 121
!!!!!!!!!!!!!!!!! end of time step  !!!!!!!!!!!!!!!!!!!!!!!!
 
 21   CONTINUE   !!!!!!!!!!!!!!!!!!!  end of calculating time
! . . . writing time at the end of program
      CALL DAT(ut,dayt,godt,verno,ut0,ut1,day,god)
      readfl = .FALSE.
      CALL ZAIT(readfl,5,nzapt,nzaps,nadrt,nadrs,isp,ldor,kdf,god,day,ut0,ut1,
     &          pole)
      CALL ZAIT(readfl,6,nzapt,nzaps,nadrt,nadrs,isp,ldor,kdf,god,day,ut0,ut1,
     &          pole)
! . . . writing time and date into inf record
      READ(4,REC=1)pole
      pole(2) = god
      pole(3) = day
      pole(4) = ut0
      pole(5) = ut1
      WRITE(4,REC=1)pole
! . . .
!  26 continue
      PRINT *, 'GSMTIP: cyclt - end'
      DEALLOCATE(vdr,potef,gins)
      END SUBROUTINE CYCLT_BAS
