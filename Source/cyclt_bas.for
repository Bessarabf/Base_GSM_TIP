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
      ! god – год, целое, формат хххх 
      ! day – день года (от 1-го января)
      ! ut0 – время, сек начала расчета
      ! utk – время, сек, окончания расчета
      ! dtt – шаг по времени, сек, в интегрировании уравнений «трубки»
      ! dts - шаг по времени, сек, в интегрировании уравнений «шара»
      ! tau - шаг по времени, сек. (предполагалось, что tau, dtt,dts могут отличаться. Сейчас все шаги равны).
      ! b - 
      ! c - 
      ! sole – массив размерности nse (15 в базовом варианте) – потоков ионизирующего излучения (КУФ) в фотонах (*1.0e11)
      ! solu – массив размерности 2*nsu (2*12) – потоков диссоциирующего излучения Солнца в фотонах (nsu) и эргах (nsu).
      ! bmpz – bz-компонента межпланентного магнитного поля (ММП)
      ! bmpy – by-компонента ММП
      ! bmod – модуль ММП,
      ! vsol – скорость солнечного ветра,
      ! csol - ,
      ! solen - массив размерности nse – потоки рассеянного КУФ
      ! solet - массив размерности nse – не используется,
      ! fa0 – индекс F10.7 текущий (среднесуточный),
      ! fs – индекс F10.7 фоновый (для модели Нусинова) или за 81 день для EUVAC. 
      ! ap0 – Ap индекс,
      ! pkp0 – Kp индекс,
      ! dst0 – Dst индекс,
      ! ae0 – AE -индекс,
      ! al0 – Al индекс,
      ! au0 – AU индекс, 
      ! nh – число точек по вертикали в «шаре» (30),
      ! ddolgs – шаг по долготе в «шаре», градус
      ! dtets – шаг по широте в «шаре», градус
      ! ddolgt – шаг по долготе в «трубке», град,
      ! ntsl - число узлов в силовой трубке
      ! nl - число силовых линий в полушарии одной долготной плоскости (17 для 5 град сетке по широте),
      ! nl2 - nl+nl+3 ,
      ! q - массив q(nv,nl)-координат дипольной системы координат (номер силовой линии, узлы вдоль силовой линии)
      ! u - массив u(nl) - координат дипольной системы координат перпендикулярной силовой линии,
      ! dtett -  – шаг по широте в «трубке», град,
      ! int - число интерполяционных параметров в “трубке” из “шара”(равно ,
      ! ins - число интерполяционных параметров в “шаре” из “трубки” (равно 6),
      ! rads - массив значений высот узлов неравномерной сетки “шара”, число элементов nh=30 , от 80 до  526 км (в сантиметрах)
      ! ntr - число узлов вдоль силовой линии дипольной системы координат,
      ! gkoor - массив размерности (2, its,ids) соответствия узлов геомагнитных координат (коширота и долгота, град) в географические.
      ! kpart - число рассчитываемых параметров в “трубке” ,
      ! kpars - число рассчитываемых параметров в “шаре” (19),
      ! nadrt - номер начальной записи параметров “трубки” в f4 ,
      ! nadrs - номер начальной записи параметров “шара” в f4 ,
      ! ns - ,
      ! idt - число долготных узлов в “трубке” (24),
      ! ids - число долготных узлов в “шаре” (24),
      ! its - число широтных узлов в шаре (37) ,
      ! nzapt - ,
      ! nzaps - ,
      ! mas - массив управляющих параметров модели (10),
      ! mast - массив управляющих параметров для вариантов расчета “трубки” (40),
      ! mass - массив управляющих параметров для вариантов расчета “шара” (30),
      ! uprt,
      ! kdf - массив номеров записей в f4, соответствующих началу параметров “шара”, потенциала, “трубки” и пр.  
      ! kdu - массив длины записей в f4,  соответствующих параметрам “шара”, потенциала, “трубки” и пр.  
      ! ldor - длина записи в f4 в байтах (4096),
      ! isp,
      ! rmaxt,
      ! par1 - массив параметров шара в одной долготной плоскости (возможно, рудимент),
      ! PAR2 - массив параметров шара в одной долготной плоскости,
      ! par - массив параметров шара в одной долготной плоскости,
      ! pari - массив размерности ni - параметры трубки в одной долготной плоскости,
      ! pole - массив размерностью 1024 для чтения параметров из файла f4,
      ! nr - размер массива параметров шара (одна долготная плоскость) =  nh*its*kpars ,
      ! verno - вспомогательная логическая переменная (.TRUE. или .FALSE.),
      ! ni - размер массива трубки ( одна долготная плоскость)  равная ntpl*int ,
      ! park,
      ! ks - число узлов в долготной плоскости “трубки” ntpl*2,
      ! nv,
      ! nsu - число спектральных интервалов ДУФ (12),
      ! nse - число спектральных интервалов КУФ (15 в базовом варианте),
      ! tet1 - коширота максимума продольных токов 1,
      ! tet2 - коширота максимума продольных токов зоны 2,
      ! tet3 - коширота максимума продольных токов зоны 3 (если есть)
      ! pdpc - разность потенциалов через шапку в кВ,
      ! eps0,
      ! om,
      ! fac1 - величина токов зоны 1,
      ! fac2 - величина токов зоны 2,
      ! fac3 - величина токов зоны 3,
      ! pglo - 4-х мерный массив, параметров “шара”, размерности (kpars,nh,its,ids)
      ! imja,
      ! dut - долгота точки, выбранной для записи параметров (вспомогательная информация),
      ! hut - высота точки  выбранной для записи параметров (вспомогательная информация),
      ! iput ,
      ! keut - длина записи в файл f4 ,
      ! koob - номер объекта (?), выбранной для записи,
      ! kut - 
      ! nara -имя входного файла с координатами объекта ,
      ! naou - имя выходного файла с координатами объекта,
      ! sut - широта объекта (?),
      ! tut - время объекта,
      ! izap ,
      ! pril - массив параметров для граничных условий “шара” (берутся из других моделей) , размерность kpa,its,ids,nt
      ! kpa - число параметров граничных условий ,
      ! nt - число временных записей для pril, обычно 12 или 24
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
