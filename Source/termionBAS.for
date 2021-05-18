! subroutine readinf, rpco, camu, wwod_D, flosu, flosuEUVAC, EUVAC, 
!     flosuN, forminf, calcdat, setka_kl, ggmraw, nachus, zerot, wwt, rtime, 
!     kodir, dat, botread
! function j52
C
c    12.12.18 Base vaiant GSM TIP based on EAGLE 
C
      program termionBAS
!      USE mo_ham_gsm

! объявление статических массивов
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Parameter(isat=8700,ldor0=4096)
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character nara(100)*25,naou(100)*25
      character imja*80

      integer uprt,uprs,verno/0/,god,day
      
      dimension mas(10),mast(40),mass(30),kdf(20),kdu(20)

!    	
      allocatable
     *            q(:),gkoor(:),pole(:)
     *           ,park(:),pari(:),par(:),par1(:),par2(:)
     *           ,pglo(:)
     *           ,pril(:,:,:,:) ! for other lowboundary model different numbers
     *           ,sole(:),solu(:),solen(:),solet(:)
     *           ,ntsl(:),u(:),rads(:)
     *           ,dut(:,:),hut(:,:),iput(:),keut(:),kut(:), 
     *            sut(:,:),tut(:,:),izap(:),nzapjet(:),lzapjet(:)

      allocate (pole(ldor0/4))
      allocate (dut(isat,100),hut(isat,100),iput(100), 
     *          keut(100),kut(100),sut(isat,100),tut(isat,100),
     *          izap(100),nzapjet(5),lzapjet(5))

      nzapjet=0
      lzapjet=(/20,20,12,12,24/) 

! объявление файлов (файла?) F4
      open(4,file='f4',status='old',access='direct',
     *     form='unformatted',recl=ldor0) ! ldor0=4096
      open(5,file='f5',status='old',access='direct',
     *     form='unformatted',recl=4096)
! чтение и разбор информационной записи
      read(4,rec=1) pole
      ldor=pole(1)
      if(ldor.ne.ldor0) then
          print*, 'incorrect length record in F4 =',ldor
          stop
	end if
	nl=pole(61)
	nl2=nl+nl+3
	allocate (ntsl(nl)) 
	call readinf(pole,ldor,god,day,ut0,ut1,utk,dtt,dts,tau,kdf,
     *                   kdu,kpart,kpars,int,ins,its,ids,idt,
     *                   ddolgt,dtett,ddolgs,dtets,ntr,nv,nh,nl,
     *                   ntsl,ntpl,ns,fa0,fs,nsu,nse,
     *                   bmpy,bmpz,vsol,csol,ap0,pkp0,dst0,
     *                   ae0,al0,au0,gamma,rmaxt,b,c,dh,rmin,
     *                   tet1,tet2,tet3,eps0,om,
     *                   fac1,fac2,fac3,pdpc,
     *                   nzapt,nzaps,nadrt,nadrs)

! выделение памяти под массивы solu и т.д.    
      allocate(sole(nse),solu(nsu),solen(nse),solet(nse))
	
! чтение danmodel
      call wwod_D(mes,god,day,ut0,utk,dtt,dts,tau,sole,solu,solet,
     *          bmpz,bmpy,bmod,vsol,csol,fa0,fs,ap0,pkp0,dst0,ae0,
     *          al0,au0,solen,nsu,nse,verno,gamma,
     *          ddolgt,dtett,nh,dh,rmin,b,c,ddolgs,dtets,
     *          dlam,rmaxt,nsill,dteta,ns,uprt,uprs,
     *          mas,mast,mass,tet1,tet2,tet3,pdpc,
     *          eps0,om,fac1,fac2,fac3,imja,kpa,ntIME,its,ids)

! сравнение на непротиворечие входных данных из информационной записи и danmodel
      if(verno.NE.0) then
              print*,' WWOD__D: incorrect input data between',
     *        'danmodel and INF record! ERRCODE=',verno
              stop
      end if
! 
! выделение памяти под динамические массивы
      ni=ntpl*int ! размер массива трубки (плоскость)
      ks=ntpl*2
      nr=nh*its*kpars	! размер массива шара (плоскость)
   
      allocate(pglo(kpars*nh*its*ids),
     *         par(nr),par1(nr),par2(nr),park(ks),pari(ni),
     *         q(nv*nl),gkoor(2*its*ids),
     *         pril(kpa,its,ids,ntIME),
     *	       rads(nh),u(nl))


! расчетная часть
!  чтение файла fiza
         if(mast(31).ne.0) then 
           call rpco(dut,hut,iput,keut,koob,kut,nara,naou,sut,tut,
     *             isat)
           open(9,file='fiza',status='old')
           read(9,'(10i5)')izap
           close(9)
         end if
  !    nse=mas(4)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    nsu=mas(5)*2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 02.03.2018
      call setka_kl(rmin,dh,gamma,ntr,nh,b,c,dtett,
     *           rmaxt,park,ks,gkoor,nv,
     *           rads,ntsl,nl,q,u,its,ids,dtets,ddolgs)

c     . . . Formed information record in FILE4
      if(uprt.eq.2) then
         call forminf(ldor,god,day,ut0,ut1,utk,dtt,dts,tau,kdf,
     *                kdu,kpart,kpars,int,ins,ddolgt,dtett,ddolgs,
     *                dtets,ntr,nv,nh,nl,ntsl,rads,mass,mast,mas,
     *                ns,fa0,fs,nsu,nse,sole,solet,solen,solu,bmpy,
     *                bmpz,vsol,csol,ap0,pkp0,dst0,ae0,al0,au0,gamma,
     *                rmaxt,b,c,dh,rmin,tet1,tet2,tet3,eps0,om,
     *                fac1,fac2,fac3,pdpc,pril,
     *                nzapt,nzaps,nadrt,nadrs)
      
      endif
      write(*,*)'**************check forminf****************',
     *ldor,god,day,ut0,ut1,utk,dtt,dts,tau,kdf	
!
c     . . . Creating global initial distribution of the parameters:
      call nachus(uprt,uprs,ut0,nh,dtets,dtett,ddolgs,ddolgt,day,
     *            ap0,fa0,fs,ntsl,nl,god,ldor,nzapt,nzaps,nadrt,
     *            nadrs,kdf,kpars,kpart,rads,park,ks,its,
     *            kdu,isp,par,pole,nr,verno,pkp0,mass,imja,mast)


	  if(verno.ne.0) then
         print *,'  nachus: incorrect input data. STOP! '
         stop
      end if


!     чтение данных приливов из фалов других моделей
      if(mass(18).eq.2) then                                 
	  call botread(day,pril,its,ids,kpa,ntIME)  
      end if  

!     . . .  Calculation cycle
      
      call cyclt_bas(god,day,ut0,utk,dtt,dts,tau,b,c,sole,solu,
     *        bmpz,bmpy,bmod,vsol,csol,solen,solet,fa0,fs,ap0,pkp0,
     *        dst0,ae0,al0,au0,nh,ddolgs,dtets,ddolgt,ntsl,nl,nl2,q,u,
     *        dtett,int,ins,rads,ntr,gkoor,kpart,kpars,nadrt,
     *        nadrs,ns,idt,ids,its,nzapt,nzaps,mas,mast,mass,
     *        uprt,kdf,kdu,ldor,isp,rmaxt,par1,PAR2,par,pari,pole,nr,
     *        verno,ni,park,ks,nv,nsu,nse,tet1,tet2,tet3,pdpc,
     *        eps0,om,fac1,fac2,fac3,pglo,imja,
     *        dut,hut,iput,keut,koob,kut,nara,naou,sut,tut,izap,
     *        pril,kpa,ntime,nzapjet,lzapjet,isat)

    
!
! освобождение памяти
      deallocate (pglo,par,par1,park,pari,PAR2,q,gkoor,pril,
     *            sole,solu,solen,solet,rads,u,ntsl,
     *            dut,hut,iput,keut,kut,sut,tut,izap,
     *            nzapjet,lzapjet,pole)  

      end program termionBAS
!------------------------------------------------------------------------------
c . . . ver. 16.04.14 - add idt to interface parameters
      subroutine readinf(pole,ldor,god,day,ut0,ut1,utk,dtt,dts,tau,kdf,
     *                   kdu,kpart,kpars,int,ins,its,ids,idt,
     *                   ddolgt,dtett,ddolgs,dtets,ntr,nv,nh,nl,
     *                   ntsl,ntpl,ns,fa0,fs,nsu,nse,
     *                   bmpy,bmpz,vsol,csol,ap0,pkp0,dst0,
     *                   ae0,al0,au0,gamma,rmaxt,b,c,dh,rmin,
     *                   tet1,tet2,tet3,eps0,om,
     *                   fac1,fac2,fac3,pdpc,
     *                   nzapt,nzaps,nadrt,nadrs)
      dimension kdf(20),kdu(20),ntsl(*),pole(1024)
      integer god,day
c        reading inform area
c
      ldor=  pole(1) 
      god =  pole(2) 
      day =  pole(3) 
      ut0 =  pole(4) 
      ut1 =  pole(5) 
      utk =  pole(6) 
      dtt =  pole(7) 
      dts =  pole(8) 
      tau =  pole(9) 
      do i=1,20
        kdf(i)=pole(i+9) 
        kdu(i)= pole(i+29)
      end do

!     число записей итерационных области шара и трубки 
!     должны совпадать с соответстующими областями шара и трубки
	if (kdu(11).ne.kdu(5)) then
	    kdu(11)=kdu(5)
	    kdu(12)=kdu(6)
	    kdu(13)=kdu(6)
	    kdu(14)=kdu(6)
	    do kk=12,15
	      kdf(kk)=kdu(kk-1)+kdf(kk-1)
	    end do
	end if
      kpart =pole(50) 
      kpars =pole(51) 
      int   =pole(52) 
      ins   =pole(53) 
      ddolgt=pole(54) 
      dtett =pole(55) 
      ddolgs=pole(56)
       
      dtets =pole(57)
      its=180.1/dtets+1
      ids=360.1/ddolgs 
	idt=360.1/ddolgt 
      ntr   =pole(58) 
      nv    =pole(59) 
      nh    =pole(60) 
      nl    =pole(61) 
      ntpl=0
      do i=1,nl
        ntsl(i)=pole(i+61)
        ntpl=ntpl+ntsl(i)
      end do
      
      ns = pole(226)
      fa0= pole(227)
      fs = pole(228)
      nsu= pole(229)
      nse= pole(230)
      
      gamma= pole(385)
      rmaxt= pole(386)
      b    = pole(387)
      c    = pole(388)
      dh   = pole(389)*1.e5
      rmin = pole(390)
      tet1 = pole(391)
      tet2 = pole(392)
      tet3 = pole(393)
      eps0 = pole(394)
      om   = pole(395)
      fac1 = pole(396)
      fac2 = pole(397)
      fac3 = pole(398)
      pdpc = pole(399)

      nzaps=  pole(600)
      nadrs=  pole(601)
      nzapt=  pole(602)
      nadrt=  pole(603)
      return
      end
!------------------------------------------------------------------------------
      subroutine rpco(dut,hut,iput,keut,koob,kut,nara,naou,sut,tut,
     *                isat)
c   out: dut  - longitude of dot orbit
c        hut  - height of dot orbit
c        iput - number of dot height profile
c        keut - length of record
c        koob - number of object
c        kut  - number of dot orbit
c        nara - name of input file object
c        naou - name of output file object
c        sut  - latitude of dot orbit
c        tut  - UT of dot orbit

       dimension dut(isat,100),hut(isat,100),iput(100),keut(100),
     *          kut(100),sut(isat,100),tut(isat,100)
       character naob(100)*4,nara(100)*25,naou(100)*25
       open(9,file='fina',status='old')
       i=0
    1 continue
         i=i+1
         read(9,'(a)',end=2) naob(i)
         goto 1
    2 continue
!      koob=i-1    ! Для NDP
       koob=i-2    ! Для FL32
c     print *, kut
       do i=1,koob
         nara(i)='raf/raf'//naob(i)
         naou(i)='ouf'//naob(i)
c       print*,naob(i),nara(i),naou(i),koob
       end do
       close(9)
       do i=1,koob
         open(10,file=nara(i),status='old')
         read(10,'(i4)') kout
c       print*,nara(i),i,kout
         kut(i)=kout
         ha=0.
         do k=1,kout
           read(10,'(i4,f9.3,2f8.2,f10.2)') ko,tut(k,i),sut(k,i),
     *                                     dut(k,i), h
           if(h.gt.ha) ha=h
           hut(k,i)=h
         end do
         close(10)
         call camu(ha,ip,kec)
         iput(i)=ip
         keut(i)=kec
c       print*,ha,ip,kec
       end do
       do i=1,koob
         kout=kut(i)
         a=tut(1,i)
         do k=2,kout
           b=tut(k,i)
           if(b.lt.a)then
             tut(k,i)=b+24.
             a=tut(k,i)
           else
             a=b
           end if
         end do
       end do
       return
       end
!------------------------------------------------------------------------------
       subroutine camu(ha,ip,kec)
       dimension he(70),hei(70)
       data hei/ 175.3, 187.8, 201.6, 216.8, 233.5, 251.8, 272.0, 294.2,
     *          318.6, 345.5, 375.0, 407.5, 443.3, 482.6, 525.9, 573.5,
     *          625.8, 683.4, 746.8, 816.4, 893.1, 977.4,1070.1,1172.1,
     *         1284.3,1407.8,1543.6,1692.9,1857.2,2037.9,2236.7,2455.4,
     *         2695.9,2960.5,3251.6,3571.7,3923.9,4311.3,4737.4,5206.2,
     *         5721.8,6289.0,6912.9,7599.1,8354.1,9184.5,10097.9,
     *        11102.7,12208.0,13423.8,14761.1,16232.3,17850.5,19630.5,
     *        21588.6,23742.4,26111.7,28717.9,31584.6,34738.1,38206.9,
     *        42022.6,46219.9,50836.9,55915.5,61502.1,67647.3,74407.0,
     *        81842.8,88557.2/
       n=70
       do i=1,n
         he(i)=hei(i)
       end do
       ip=j52(n,ha,he)+2
       if(ip.gt.70)ip=70
       kec=(ip*8+480+9)*4
       return
       end
!------------------------------------------------------------------------------
       function j52(n,u,x)
       dimension x(n)
       i=n
       s=x(i)
       i=n-1
       if(u.ge.s) goto 2
       i=1
       s=x(i+1)
       if(u.lt.s) goto 2
       j=n+1
    1 continue
         k=(i+j)/2
         s=x(k)
         if(u.lt.s) j=k
         if(u.ge.s) i=k
       if(j.gt.i+1) goto 1
    2 continue
       j52=i
       return
       end
!------------------------------------------------------------------------------
c    ver     21.12.2012 - dtets,ddolgs & dtett,ddolgt - read from danmodel
c            and then compare with inf record  
c            25/05/2012 - pril with kpa & ntIME
c    version 25.03.2012 bmod - modulus of solar magnetic field
c    version 11.03.2010 flosu - flosuN
c    
      subroutine wwod_D(mes,god,day,ut0,utk,dtt,dts,tau,sole,solu,solet,
     *          bmpz,bmpy,bmod,vsol,csol,fa0,fs,ap0,pkp0,dst0,ae0,
     *          al0,au0,solen,nsu,nse,verno,gamma,
     *          ddolgt,dtett,nh,dh,rmin,b,c,ddolgs,dtets,
     *          dlam,rmaxt,nsill,dteta,ns,uprt,uprs,
     *          mas,mast,mass,tet1,tet2,tet3,pdpc,
     *          eps0,om,fac1,fac2,fac3,imja,kpa,ntIME,its,ids)
      dimension sole(nse),solu(nsu),solen(nse),solet(nse)
      dimension mas(10),mast(40),mass(30)

      dimension dlw1(20),dlw2(20),dw1(100),dw2(100)
      integer uprt,uprs,god,day,verno,uth,utm,tkh,tkm
      integer mmes(12)/31,28,31,30,31,30,31,31,30,31,30,31/,b1,c1
      character *8 mes1(12),imja*80
      open(7,file='danmodel',status='old')
c    *    ,form='formatted',recl=80)
      data mes1/'january ','february','march   ','april   ',
     *     'may     ','june    ','july    ','august  ','septemb.',
     *     'oktober ','november','december'/
c
  902 format(' ***************************************',
     *       '***************************************'/
     *       ' *',20x,'MODEL  GSM-TIP',42x,'*'/
     *       ' ***************************************',
     *       '***************************************')
  300 format(//a80)
  901 format(//2x,i4,2x,i2,2x,i2,4x,i2,2x,i2,3x,i2,2x,i2,
     *       4x,f5.2,5x,f5.2,3x,f5.2)
  910 format('   god  = ',i4,3x,
     *       '   day  = ',i2,1x,a8,'  - ','day of god(',i3,')'/
     *       '   ut0  = ',i2,' hour  ',i2,' min.',2x,'(',f8.0,' sek.)'/
     *       '   tk   = ',i2,' hour  ',i2,' min.',2x,'(',f8.0,' sek.)'/
     *       '   tau  = ',f5.1,'(min.) ',2x,'(',f5.0,'sek. )'/
     *       '   dts  = ',f5.1,'(min.)',2x,' (',f5.0,' sek.)'/
     *       '   dts  = ',f5.1,'(min.)',2x,' (',f5.0,' sek.)')    
      open(10,file='fort.10')
      write(10,902)
      print 902

      read (7,300)imja
      read(7,901)god,day,mes,uth,utm,tkh,tkm,tau,dts,dtt
      taum=tau
      tau=tau*60
      dtsm=dts
      dts=dts*60
      dttm=dtt
      dtt=dtt*60
      nc1=god/4
      nc1=god-nc1*4
      if(nc1.eq.0)mmes(2)=29
      nden=day
      if(mes.ne.1) then ! day of year
         nc=mes-1
         do i=1,nc
            day=day+mmes(i)
         end do
      end if 
      ut0=60*(uth*60+utm)
      tk=60*(tkh*60+tkm)
      k=mes
      print 910,god,nden,mes1(k),day,uth,utm,ut0,
     *       tkh,tkm,tk,taum,tau,dtsm,dts,dttm,dtt
c      write(10,910)god,nden,mes1(k),day,uth,utm,ut0,
c     *       tkh,tkm,tk,taum,tau,dtsm,dts,dttm,dtt

      read(7,923)ns,uprs,uprt
c     write(10,944)ns,uprt,uprs
  944 format('   ns   = ',i2,'   uprt = ',i1,'   uprs = ',i1)
  923 format(/3x,i2,5x,i1,4x,i1)

      read(7,921)mass

      read(7,1921)mast

      read(7,921)mas

      write(10,922)mass,mast,mas
  921 format(/1x,i3,19i4/1x,i3,9i4)
 1921 format(/1x,i3,19i4/1x,i3,19i4) 
  922 format(15x,' MASS :'/20i4/10i4/
     *       15x,' MAST :'/20i4/20i4/
     *       15x,' MAS :'/10i5)
      if(mas(5).gt.20)go to 120
      if(mas(4).gt.100)go to 120

      utk=ut0+tk
!  fa0 - background F10.7, fs - daily F10.7
      read(7,919)fa0,fs
	   
      nn=mas(4)
	
      read(7,903)(dw1(i),dw2(i),sole(i),i=1,nn)
  919 format(/1x,2(f4.0,1x)) 
  903 format(/3(f7.1,1x,f7.1,1x,f10.5,1x))
      
      if(mas(8).ne.0) then
          if (mas(8).eq.1 ) then
     	     call flosu(god,fa0,sole,nn)
	     print*,' flosu'
	    else if (mas(8).eq.2) then
     	      call flosuN(fs,fa0,sole,nn)
	      print*,' flosuN'
          else if (mas(8).eq.3) then
     	      call flosuEUVAC(fs,fa0,sole,nn)
          else
            print*,' wwod: incorrect mas(8)'
            stop
 	    end if
      end if
	
      m=mas(4)/2
      mr=mas(4)-m*2
      nc1=m
      nc2=m
      nc3=m
      
      if(mr.ge.1)nc1=m+1
      if(mr.ge.2)nc2=m+1
      if(mr.ge.3)nc3=m+1
    
      nc2=nc2+nc1
      nc3=nc2+nc3
      write(10,932)fa0,fs
      write(10,933)
  932 format(/'   fa0 =   ',f4.0,'    fs  =   ',f4.0)
  933 format('             SOLE '/
     *      (2('     длина волны ','       sole  ')))
  934 format('  ',2(f7.1,'-',f7.1,2x,f10.5,3x))
      if(mas(4).lt.2) go to 21
      do k=1,m
         write(10,934)dw1(k),dw2(k),sole(k),
     *                dw1(k+nc1),dw2(k+nc1),sole(k+nc1)
      end do
   21 continue
      if(mr.eq.0) go to 8
      m=m+1
      mr=mr*m
      write(10,934)(dw1(k),dw2(k),sole(k),k=m,mr,m)
    8 continue
      go to 130
  120 continue
      write(10,800)mas(4),mas(5)
      print 800,mas(4),mas(5)
  800 format(' ********   mas(4)=',i7,
     *       '  mas(5)=',i7,' ********')
      stop
  130 continue
     
      read(7,929)(solet(i),i=1,nn)
      write(10,777)
      write(10,976)(solet(i),i=1,nn)
  777 format('             SOLET')
  929 format(/3(3x,7f8.3/))
  976 format(' ',2(7(f6.3,3x)/' '))
     
      m=mas(4)/2
      mr=mas(4)-m*2
      read(7,927)(solen(i),i=1,nn)
      write(10,771)
  927 format(/2(7(f10.5,1x)/),f10.5)
  771 format(10x,'   SOLEN'/
     *       2('     длина волны ','       solen '))
   
      if(mas(4).lt.2) go to 22
         write(10,*)  solen
   22 continue
C     if(mr.eq.0) go to 88
C     m=m+1
C     mr=mr*m
C     write(10,934)(dw1(k),dw2(k),solen(k),k=m,mr,m)
C  88 continue
      nn=mas(5)
      read(7,917)(dlw1(i),dlw2(i),solu(i),solu(i+nn),i=1,nn)
  917 format(/2(2(f7.1,1x),2(f10.5,1x)))
      m=mas(5)
      mr=mas(5)-m
      nc1=m
      nc2=m
      if(mr.eq.0)go to 6
      if(mr.ge.1)nc1=m+1
    6 continue
      write(10,940)(dlw1(k),dlw2(k),solu(k),solu(k+nn),
     *       k=1,m)

      read(7,905)bmpy,bmpz,bmod,vsol,csol
      write(10,920)bmpy,bmpz,bmod,vsol,csol
      read(7,907)ap0,pkp0,dst0,ae0,al0,au0
      write(10,924)ap0,pkp0,dst0,ae0,al0,au0

  940 format(10x,'  SOLU  '/'     длина волны ',
     * 9x,'solu: fot/(cm**2*sek) ','     solu: erg/(sm**2*sek)  '/
     * 20(/'   ',f7.1,'-',f7.1,10x,f10.5,15x,f10.5,7x))

  905 format(/5(1x,f7.2))
  920 format('   bmpy =   ',f7.2,1x,'   bmpz =   ',f7.2,1x,
     *       '   bmod=   ',f7.2,1x,
     *       '   vsol =', f 7.2,1x,'   csol =  ',f7.2)
  907 format(/3(1x,f5.0),7x,3(1x,f5.0))
  924 format('   ap0 = ',f6.0,5x,' pkp0 = ',f6.0,3x,'  dst0 =   ',f6.0/
     *       '   ae0 = ',f6.0,5x,' al0  = ',f6.0,3x,'  au0  =   ',f6.0)

! чтение и проверка на совпадение в инф.записи шагов трубки 
      read(7,909)ntr,gamma,ddolgt0,dtett0,rmaxt,b1,c1
      if(ddolgt0.ne.ddolgt) verno=1
      if(dtett0.ne.dtett) verno=2
  909 format(/1x,i3,1x,f8.4,3x,f4.0,2x,f4.0,4x,f5.2,2x,i3,1x,i3)
      
      write(10,926) gamma,rmaxt,b1,c1
  926 format(/' gamma=',f4.1,'; rmaxt=',f5.2,'*re;',2x,
     *        ' b=',i3,'; c=',i3)

      b=b1
      c=c1
      re=6371.02
      rmaxt=rmaxt*re

!  чтение и проверка на совпадение шагов шара и трубки
      read(7,925)nh,dh,rmin,ddolgs0,dtets0,dlam,nsill,dteta
      if(floor(ddolgt0).ne.floor(ddolgt)) verno=3
      if(floor(dtett0).ne. floor(dtett)) verno=4 
      write(10,964)nh,dh,rmin,ddolgs,dtets,dlam,nsill,dteta
 925  format(/1x,i3,1x,f4.0,3(3x,f4.0),11x,f4.0,6x,i2,3x,f4.0)
 964  format(/'  nh     = ',i3,'    dh    = ',f4.0,'  rmin  =  ',f4.0,/
     *        '  ddolgs = ',f4.0,'   dtets = ',f4.0//
     *        '  dlam   = ',f4.0,'   nsill = ',i4,'  dteta = ',f4.0)
 
      read (7,931)tet1,tet2,tet3,eps0,om,fac1,fac2,fac3,pdpc
  931 format(/2x,f3.0,3x,f3.0,3x,f3.0,2x,e6.2,2x,f4.2,4(2x,e6.2))
      write(10,939)tet1,tet2,tet3,eps0,om,fac1,fac2,fac3
  939 format(/'   tet1 = ',f4.0,3x,'tet2 = ',f4.0,5x,'tet3 = ',f4.0/
     *        '   eps0 = ',1pe9.2,7x,'om = ',0pf5.2/
     *        '   fac1 = ',1pe9.2,3x,'fac2 = ',1pe9.2,
     *        '   fac3 = ',1pe9.2)

!  чтение параметров для задания массива PRIL
	read(7,*) 
	

      read(7,*) kpa,ntime

      close(10)
      close(7)
      return
      end
!------------------------------------------------------------------------------
      subroutine flosu(godi,fa0,sole,mas4)
      dimension al(35),bl(35),sole(mas4),bla(3),blr(3),
     *          sol(38),soler(8),ale(8),d(8),solei(8)
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
      god=godi
c
      gg=(god-gb)/t
      if(gg.gt.3.) gg=gg-3.
      if(gg.gt.2.) gg=gg-2.
      if(gg.gt.1.) gg=gg-1.
      ss=abs(sin(pi*gg))
c . . . F10.7 - минимальный = 63
      ff=63.
      if(ss.ge.1.e-3)ff=63.+482.*(ss)**(+3.7)*exp(-5.2*gg)
c
      print*,'flosu ff=',ff
      fff=(ff-60.)**2
      fff=fff**0.333
      ffa=(fa0-ff)**2
      ffa=ffa**0.333
      fir=blr(1)+blr(2)*fff+blr(3)*ffa
c
      do 1 i=1,35
        sol(i+2)=(al(i)+bl(i)*fir)*fir
    1 continue
      sol(38)=bla(1)+bla(2)*fff+bla(3)*ffa
c
      ri820=0.023*fa0-1.440
      do 2 i=1,8
        d(i)=1.56/ale(i)+0.22
        solei(i)=soler(i)*(ri820/soler(1))**d(i)*1.e-2
    2 continue
      sol(1)=solei(1)+solei(2)+solei(3)
      sol(2)=0.
      do 3 i=4,8
        sol(2)=sol(2)+solei(i)
    3 continue
      sol(2)=sol(2)
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
      sole(15)=sol(38)
      return
      end
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
      subroutine flosuN(ff,fa0,sole,mas4)
      dimension al(35),bl(35),sole(mas4),bla(3),blr(3),
     *          sol(38),soler(8),ale(8),d(8),solei(8)
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
      if(ff.lt.63) ff=63.
c
      fff=(ff-60.)**2
      fff=fff**0.333
      ffa=(fa0-ff)**2
      ffa=ffa**0.333
      fir=blr(1)+blr(2)*fff+blr(3)*ffa
c
      do 1 i=1,35
        sol(i+2)=(al(i)+bl(i)*fir)*fir
    1 continue
      sol(38)=bla(1)+bla(2)*fff+bla(3)*ffa
c
      ri820=0.023*fa0-1.440
      do 2 i=1,8
        d(i)=1.56/ale(i)+0.22
        solei(i)=soler(i)*(ri820/soler(1))**d(i)*1.e-2
    2 continue
      sol(1)=solei(1)+solei(2)+solei(3)
      sol(2)=0.
      do 3 i=4,8
        sol(2)=sol(2)+solei(i)
    3 continue
      sol(2)=sol(2)
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
      sole(15)=sol(38)
      return
      end
!------------------------------------------------------------------------------
      subroutine forminf(ldor,god,day,ut0,ut1,utk,dtt,dts,tau,kdf,
     *                   kdu,kpart,kpars,int,ins,
     *                   ddolgt,dtett,ddolgs,dtets,ntr,nv,nh,nl,
     *                   ntsl,rads,mass,mast,mas,ns,fa0,fs,nsu,nse,
     *                   sole,solet,solen,solu,
     *                   bmpy,bmpz,vsol,csol,ap0,pkp0,dst0,
     *                   ae0,al0,au0,gamma,rmaxt,b,c,dh,rmin,
     *                   tet1,tet2,tet3,eps0,om,
     *                   fac1,fac2,fac3,pdpc,pril,
     *                   nzapt,nzaps,nadrt,nadrs)
      dimension kdf(20),kdu(20),ntsl(44),pole(1024),pril(*)
      dimension sole(40),solet(40),solen(40),solu(24)
      dimension mas(10),mass(30),mast(40),rads(30)
      integer god,day
      character*8 namemes,im*80
  900 format(a25)
  901 format(3(5i6/),4(10i6/),2i6,2e10.3)
c        Формирование информационной области
c
      pole(1)  = ldor
      pole(2)  = god
      pole(3)  = day
      pole(4)  = ut0
      ut1=ut0
      pole(5)  = ut1
      pole(6)  = utk
      pole(7)  = dtt
      pole(8)  = dts
      pole(9)  = tau
      do i=1,20
        pole(i+9) =kdf(i)
        pole(i+29)=kdu(i)
      end do
      pole(50) = kpart
      pole(51) = kpars
      pole(52) = int
      pole(53) = ins
      pole(54) = ddolgt
      pole(55) = dtett
      pole(56) = ddolgs
      pole(57) = dtets
      pole(58) = ntr
      pole(59) = nv
      pole(60) = nh
      pole(61) = nl
      do i=1,nl
        pole(i+61)=ntsl(i)
      end do
      do i=1,nh
        pole(i+105)=rads(i)
      end do
      do i=1,30
        pole(i+135)=mass(i)
      end do
      do i=1,40
        pole(i+165)=mast(i)
      end do
      do i=1,10
        pole(i+205)=mas(i)
      end do
      pole(226)= ns
      pole(227)= fa0
      pole(228)= fs
      pole(229)= nsu
      pole(230)= nse
      do i=1,nse
        pole(i+230) = sole(i)
        pole(i+270) = solet(i)
        pole(i+311) = solen(i)
      end do
      do i=1,nsu
        pole(i+350) = solu(i)
      end do
      pole(375)= bmpy
      pole(376)= bmpz
      pole(377)= vsol
      pole(378)= csol
      pole(379)= ap0
      pole(380)= pkp0
      pole(381)= dst0
      pole(382)= ae0
      pole(383)= al0
      pole(384)= au0
      pole(385)= gamma
      pole(386)= rmaxt
      pole(387)= b
      pole(388)= c
      pole(389)= dh*1.e-5
      pole(390)= rmin
      pole(391)= tet1
      pole(392)= tet2
      pole(393)= tet3
      pole(394)= eps0
      pole(395)= om
      pole(396)= fac1
      pole(397)= fac2
      pole(398)= fac3
      pole(399)= pdpc
      do i=1,143
        pole(i+399) = pril(i)
      end do
      sum=0
      pole(600)= nzaps
      pole(601)= nadrs
      pole(602)= nzapt
      pole(603)= nadrt
      do i=1,20
      sum = sum+kdu(i)
      end do
      pole(1000) =sum
      pole(1004) = 4
      pole(1005) = 5
      pole(1006) = 6
      pole(1007) = 7
      pole(1008) = 8
      pole(1009) = 9
      pole(1010) = 10
      pole(1011) = 11
      pole(1012) = 12
      pole(1013) = 13
      pole(1014) = 14
      write(4,rec=1)pole
      call calcdat(ut0,ut1,day,god,ng,nden,namemes,np,min)
  920 format(i19,' year',i3,1x,a8,i3,' hour.',i2,' min.'/)
      print   920,ng,nden,namemes,np,min
	print*,' god,day ', god,day

      return
      end
!------------------------------------------------------------------------------
c                Расчет даты и времени
      subroutine calcdat(ut0,ut1,day,god,ng,nden,namemes,np,min)
      integer god,day
      integer kdm(12)/31,29,31,30,31,30,31,31,30,31,30,31/
      character*8 mes(12),namemes
      data mes/'january ','february','march   ','april   ',
     *         'may     ','june    ','july    ','august  ',
     *         'septembr',
     *         'october ','november','december'/
  920 format(i19,' year',i3,1x,a8,i3,' hour.',i2,' min.'/)
c                 Расчет времени суток, дня и года
      nd=day
      ut=ut1
      ng=god
      ngd=ng/4
      ngd=ng-ngd*4
      ngod=365
      if(ngd.eq.0)ngod=366
      i=0
    1 if(ut.lt.86400)go to 2
        ut=ut-86400
        i=i+1
        go to 1
    2 continue
      nd=nd+i
      if(nd.le.ngod)go to 3
        nd=nd-ngod
        ng=ng+1
    3 continue
c                 Расчет дня месяца и номера месяца
      np=ng/4
      np1=ng-np*4
      if(np1.eq.0)go to 4
        kdm(2)=28
    4 continue
      nmes=1
      nden=nd
      if(nden.le.31)go to 7
        do 6 i=1,12
          nden=nden-kdm(i)
          nmes=nmes+1
          if(nden.gt.kdm(i+1))go to 5
            go to 7
    5     continue
    6   continue
    7 continue
c          nmes - номер месяца
c          nden - день месяца
      nsek=0
      np=0
      min=0
      nut0=ut
      np=nut0/3600
      np1=nut0-np*3600
      min=np1/60
c     nsek=np1-min*60
      namemes=mes(nmes)
c     print   920,ng,nden,namemes,np,min
      return
      end
!------------------------------------------------------------------------------
c                S E T K A
c ver 02.03.2018 
      subroutine setka_kl(rmin,dh,gamma,ntr,nh,b1,c1,dtett,
     *                 rmaxt,par,nr,gkoor,nv,
     *                 rads,ntsl,nl,q,u,its,ids,dtets,ddolgs)
      dimension par(nr),rads(nh),
     *          q(nv,nl) ,u(nl),ntsl(nl),gkoor(2,its,ids)
      double precision pi,a,b,c,d,f,g,p,re,e
      data  pi/3.14159265359d0/
      b=b1
      c=c1
      re=6.37102d8
      a=pi/180.d0
      d=-dtett
      htm=rmaxt*1.e5-re
      hm=rmin*1.e5
      dh=dh*1.e5
      rads(1)=hm
c  . . . высоты "шара"
      do i=2,nh
        rads(i)= rads(i-1)+dh*gamma**(i-2)
      end do
      ns=0
      e=rads(ntr)
      nl=0
C    . . . сетка "трубки":
C    . . . цикл по коширотам от юж.полюса до экватора
    2 continue
      if(b.lt.c)go to 14
        nl=nl+1
        li=ns+1
        par(li)=e
        f=dsin(b*a)
        f=f*f/(re+e)
        par(li+1)=b*a
        g=e+dh*gamma**(ntr-1)
        h=1.d0/f-re           ! . . . высота вершины линии
        o=amin1(h,htm)        ! . . . высота верхнего узла линии
        li=li+2
        i=1
c . . . Цикл по I с предусловием
    3   if(g.gt.o)go to  4
C   . . .  Восходящая половина линии
          i=i+1
          par(li)=g
          p=(re+g)*f
          par(li+1)=pi-datan(dsqrt(p/(1.d0-p)))
          g=g+dh*gamma**(i+ntr-2)
          li=li+2
          go to 3
    4   continue
c . . . Конец цикла по I
        j=i+1
        if(o.eq.h) then
c . . . Замкнутая линия:
          par(li)=h              ! . . . экваториальная точка
          par(li+1)=pi*.5
          nx=j+j-1
          do k=1,i
            m=ns+2*(j+k)-1
            l=ns+2*(j-k)-1
            par(m)=par(l)
            par(m+1)=pi-par(l+1)
! ...широты основания линий - одинаковы! 02.03.2018
            if(k.eq.i)par(m+1)=(180.-b)*a 
          end do
        else
c . . . Разомкнутая линия:
          par(li)=o
          g=(re+o)*f
          h=pi-datan(dsqrt(g/(1.d0-g)))
          par(li+1)=h
          li=li+2
          par(li)=o
          par(li+1)=pi-h
          nx=j+j
          do k=1,j
            l=ns+2*(j-k)+1
            m=ns+2*(j+k)-1
            par(m)=par(l)
            par(m+1)=pi-par(l+1)
! ...широты основания линий - одинаковы! 02.03.2018
            if(k.eq.j)par(m+1)=(180.-b)*a
          end do
        end if
        li=ns+1
        f= sin(par(ns+2))
        u(nl)=re/(re+par(li))*f*f
        do i=1,nx
          o=par(li)
          f=par(li+1)
          h=re/(re+o)
          q(i,nl)=h*h*dcos(f)
          li=li+2
        end do
        i=nx/2
        j=(nx+1)/2
        if(i.ne.j) q(j,nl)=0.
        ntsl(nl)=nx
        ns=ns+nx*2
        b=b+d
        go to 2
   14 continue
c . . . Преревод в градусы коширот
      do i=2,ns,2
        par(i)=par(i)/a
      end do
      nsp=ns/2
      print 900,nsp,nl,(ntsl(i),i=1,nl)
  900 format(' SETKA:   sum : ntsl(nl) =',i5,' lines = ',i4/
     *       10i5/10i5)
      do i=1,ids
        dolm=ddolgs*(i-1)
        do j=1,its
            tetm=dtets*(j-1)
          call ggmraw(1,dolg,tet,dolm,tetm)
          gkoor(1,j,i)=tet
          gkoor(2,j,i)=dolg
        end do
      end do
      return
      end
!------------------------------------------------------------------------------
      subroutine ggmraw(art,dolg,tet,dolm,tetm)
! art - var of convertation 
! 0 - geo to mag coordinate
! 1 - mag to geo coordinate 
! input and output in gratitudes
      integer art
      double precision
     *    zpi,faktor,cbg,ci,si,xlm,bm,cbm,sbm,
     *    clm,slm,sbg,bg,slg,clg,xlg,ylg
      zpi=6.28318530718
      faktor=0.01745329252
      cbg=11.4*faktor
      ci=dcos(cbg)
      si=dsin(cbg)
      if(art.eq.0) go to 10
        xlm=dolm
        bm=90.-tetm
        cbm=dcos(bm*faktor)
        sbm=dsin(bm*faktor)
        clm=dcos(xlm*faktor)
        slm=dsin(xlm*faktor)
        sbg=sbm*ci-cbm*clm*si
        bg=dasin(sbg)
        cbg=dcos(bg)
        slg=(cbm*slm)/cbg
        clg=(sbm*si+cbm*clm*ci)/cbg
        if(clg.gt.1..or.(1.-clg).lt.1.e-10)clg=1.
        if(clg.lt.-1..or.(1.+clg).lt.1.e-10)clg=-1.
        xlg=dacos(clg)
        if(slg.lt.0.0) xlg=zpi-dacos(clg)
        bg=bg/faktor
        xlg=xlg/faktor
        xlg=xlg-69.8
        if(xlg.lt.0.0) xlg=xlg+360.0
        tet=90.-bg
        dolg=xlg
        go to 20
   10   bg=90.-tet
        xlg=dolg
        ylg=xlg+69.8
        cbg=dcos(bg*faktor)
        sbg=dsin(bg*faktor)
        clg=dcos(ylg*faktor)
        slg=dsin(ylg*faktor)
        sbm=sbg*ci+cbg*clg*si
        bm=dasin(sbm)
        cbm=dcos(bm)
        slm=(cbg*slg)/cbm
        clm=(-sbg*si+cbg*clg*ci)/cbm
        if(clm.lt.-1.0) clm=-1.0 
	  if(clm.gt.1.0) clm=1.0 
	  xlm=dacos(clm)
        if(slm.lt.0.0) xlm=zpi-dacos(clm)
        bm=bm/faktor
        xlm=xlm/faktor
        dolm=xlm
        if(abs(tet-180.).lt.1.e-3)dolm=0.
        tetm=90.-bm
   20 continue
      return
      end
!------------------------------------------------------------------------------
c    . . . Вариант начала счета (определяется UPRS & UPRT)
      subroutine nachus(uprt,uprs,ut0,nh,dtets,dtett,ddolgs,ddolgt,
     *                 day,ap0,fa0,fs,ntsl,nl,god,ldor,nzapt,nzaps,
     *                 nadrt,nadrs,kdf,kpars,kpart,rads,park,ks,its,
     *                 kdu,isp,par,pole,nr,verno,pkp0,mass,imja,mast)
      dimension ntsl(nl),par(nr),pole(ldor/4),rads(nh),mast(40),
     *          kdf(20),kdu(20),park(ks),mass(30)
      character *1 imja(80)
      integer uprt,uprs,day,god,verno
      logical readfl
      mdor=ldor/4
c                         Shar
      nfw=11
      nfr=5
      nfile=5
      if(uprs.eq.0)then ! нулевые начальные условмя - не применяются
c
c        call zeros(day,ap0,fa0,fs,ut0,nh,dtets,ddolgs,rads,
c     *       ldor,kdf,isp,kpars,par,pole,nr,verno,its,pkp0,mass)
        print   921
  921   format(' zeros(shar) -> Nachus  ')
c       write (10,921)
        readfl=.false.
        call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *            god,day,ut0,pole,verno,imja)
      else if(uprs.eq.2) then ! время счета задается в danmodel 
c           write (10,941)
            print   941
  941       format(' Disk(shar) -> Nachus   ')
            readfl=.true.  ! читаем текущее время и дату в file4
            call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *                nt,nt1,unt,pole,verno,imja)
            readfl=.false. ! пишем дату и время   
            call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *                god,day,ut0,pole,verno,imja)
      else if(uprs.eq.3) then ! продолжение счета
    
            print*, ' Calculation;   uprt = uprs = 3 '
            readfl=.true.  ! читаем дату и время 
            call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *                nt,nt1,unt,pole,verno,imja)
           ! go to 10
      else 
            print*, ' Incorrect uprs -> stop !   '
            stop      
      end if
      if(verno.eq.1) return
c
c                       Trubka
      nfw=12
      nfr=6
      nfile=6
      readfl=.false.
      if(uprt.eq.0) then
            call zerot(ddolgt,ntsl,nl,ut0,kdf,ldor,isp,kpart,
     *      par,pole,nr,park,ks,mast)
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *       god,day,ut0,pole,verno,imja)
            print*,' zerot(trubka) -> Nachus '
      else if(uprt.eq.1) then 
c           call mldisk(nfile,kdf,kdu,ldor,isp,par,nr,verno)
            readfl=.true.
            print*, ' ML(trubka) ->  Nachus '
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 nt,nt1,unt,pole,verno,imja)
            readfl=.false.
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 god,day,ut0,pole,verno,imja)
      else if(uprt.eq.2) then 
            print*, ' Disk(trubka) -> Nachus '
            readfl=.true.
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 nt,nt1,unt,pole,verno,imja)
            readfl=.false.
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 god,day,ut0,pole,verno,imja)
      end if
      return
      end
!------------------------------------------------------------------------------
      subroutine zerot(ddolgt,ntsl,nl,ut0,kdf,ldor,isp,kpart,
     *                 par,pole,nr,park,ks,mast)
       integer art
      logical readfl
      dimension kdf(20),par(nr),pole(ldor/4),msum(82)
      dimension ntsl(nl),park(ks),mast(40)
      real hmi(3)/300.,700.,500./
!      real cim(3)/5.e5,1.e1,1.e-3/
      real hi(3)/40.,400.,200./
      zbaz=125.
      tzbaz=300.
      tbesk=1000.
      s=0.014
      re=6371.02
      pi=3.1415926
      cr=pi/180.
      msum(1)=0
      do 10 i=2,nl
        msum(i)=msum(i-1)+ntsl(i-1)
   10 continue
cc    nfile=3
cc    kpar=2
cc    readfl=.true.
cc    dolm=0.
cc    md=1
cc    call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
cc   *       kdf,ldor,isp,
cc   *       md,park,pole,nr,mast)
      dolm=0.
    5 if(dolm.ge.360.)go to 1
        nomsl=1
   4    if(nomsl.gt.nl)go to 2
        l=kpart*msum(nomsl)
        l1=2*msum(nomsl)
        nt=ntsl(nomsl)
        do 3 i=1,nt
          m1=l1+2*(i-1)+1
          alt=park(m1) /1.e5
          m=l+kpart*(i-1)+1
          tet=park(m1+1)
          art=1
          call ggmraw(art,dol,tetg,dolm,tet)
          d=dol*cr+pi/2.
          pr=0.55+0.45*sin(d)
          do 6 k=1,3
            cnt=(alt-hmi(k))/hi(k)
            if(cnt.gt.70.)cnt=70.
c
ccc         par(m)=cim(k)*exp(1.-cnt-exp(-cnt))*pr
ccc         if(par(m).lt.1.e-25)par(m)=1.e-25
            par(m)=1.e-3
            m=m+1
    6     continue
          do 7 k=4,6
            par(m)=0.
            m=m+1
    7     continue
c
           sss=s*(alt-zbaz)
           if(sss.gt.70.)sss=70.
          tn=tbesk-(tbesk-tzbaz)*exp(-sss)
c         tn=tbesk-(tbesk-tzbaz)*exp(-s*(alt-zbaz))
          par(m)=tn
cc        par(m)=1000.
          par(m+1)=tn
cc        par(m+1)=1000.
    3   continue
        nomsl=nomsl+1
        go to 4
    2 continue
        kpar=kpart
        readfl=.false.
        md=1
        nfile=6
        call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
     *       kdf,ldor,isp,
     *       md,par,pole,nr,mast)
c
        dolm=dolm+ddolgt
        go to 5
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine wwt (readfl,nfile1,kpar,dolm,ddolgt,nomsl,
     *       ntsl,nl,kdf,ldor,isp,md,par,pole,nr,mast)
      dimension ntsl(nl),kdf(20),par(nr),pole(ldor/4),mast(40)
      logical readfl
      nfile=nfile1
      if(nfile.eq.3)print 101,nfile,readfl
  101 format(' wwt :   nfile=',i3,'  readfl=',l8)
      if(nfile.eq.14)go to 8
      go to 9
    8 continue
      if(mast(13).eq.0) nfile=13
c     if(mast(1).eq.1)nfile=13
    9 continue
      nob1=0
      nob2=0
      ns=0
      do 1 i=1,nl
        ns=ns+ntsl(i)
    1 continue
      if(nfile.eq.3)go to 3
        if(dolm.eq.0.)go to 2
          nn=dolm/ddolgt
          nob2=kpar*ns*nn
    2   continue
    3 continue
      if(md.ne.0)go to 6
        if(nomsl.eq.1)go to 5
          ns1=0
          k=nomsl-1
          do 4 i=1,k
            ns1=ns1+ntsl(i)
    4     continue
          nob1=kpar*ns1
    5   continue
        lpar=ntsl(nomsl)*kpar
        go to 7
    6 continue
        lpar=kpar*ns
    7 continue
      nob=nob1+nob2
      call inpout(readfl,nfile,kdf,ldor,
     *       isp,nob,lpar,par,pole,nr)
      return
      end
!------------------------------------------------------------------------------
c              Чтение или запись даты и времени
c     ver 29/05/2014
 
      subroutine rtime(readfl,nfile,nzap,adres,isp,ldor,kdf,uprt,
     *                god,day,ut0,pole,verno,imja)
      logical readfl

      integer uprt,verno,god,day,adres

      character*1 ntexto(80),imja(80)
      dimension pole(ldor/4),imjan(80),imjar(80),kdf(20)

      integer nfa(20)/4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,
     *       0,0,0,0,0/

      integer kdm(12)/31,28,31,30,31,30,31,31,30,31,30,31/
      character*8 mes(12)
      data mes/'january ','february','march   ','april   ',
     *         'may     ','june    ','july    ','august  ',
     *         'septembr','october ','november','december'/
      mdor=ldor/4

      isp=kdf(nfile)+nzap

      nf=nfa(nfile)
	print*,nf,isp,nfile,kdf(nfile),nzap
	
      read(nf,rec=isp)pole
 	print*,nf,isp,pole(adres+3)
	
      if(.not.readfl) then ! write data 
         pole(adres+3)=god
         pole(adres+2)=day
         pole(adres)=ut0
         pole(adres+1)=ut0

! . . . преобразуем название месяца 
         call kodir(imja,imjan,0,imjar,ntexto)
         jjj=3
         do iii=1,80
            pole(adres+jjj+iii)=imjan(iii)
         end do
! пишем название месяца за датой
        isp=kdf(nfile)+nzap
        write(nf,rec=isp) pole
      else  ! читаем дату и время  go to 8
    
        ut0=pole(adres)
        ut1=pole(adres+1)
        day=pole(adres+2)
        god=pole(adres+3)

        jjj=3

        do iii=1,80
          imjar(iii)=pole(jjj+iii+adres)
        end do 
! . . .преобразуем в текст
        call kodir(imja,imjan,1,imjar,ntexto)
        if(uprt.ge.2)then
          ut=ut0
          nd=day
          ng=god
          call dat(ut0,day,god,verno,ut,ut1,nd,ng)
        end if
     
        if(verno.eq.1) return ! end of rtime
! . . . проверка на високосность          
          np=god/4
          np1=god-np*4
          kdm(2)=29
          if(np1.ne.0) kdm(2)=28 ! february
     
          nmes=1
          nden=day
! . . . ищем месяц и день месяца:
          if(nden.gt.31) then 
          
			do while(nden.gt.kdm(nmes)) 
				nden=nden-kdm(nmes)
				nmes=nmes+1
				if(nmes.gt.12) then 
				  print*, ' incorrect number day=',day,' STOP'
				  stop
				end if 
			end do 
	    end if
c     nmes- nomer mesjaca
c     nden- den  mesjaca
          nsek=0
          np=0
          min=0
          nut0=ut0
          np=nut0/3600
          np1=nut0-np*3600
          min=np1/60
          nsek=np1-min*60
          print   920,god,nden,mes(nmes),np,min,nsek,ntexto
    9   continue
      end if

  920 format(' *',i5,' year',i3,1x,a8,' time: ',
     *   i3,' hour.',i2,' min.',i2,' sek. ,',8a1/' ',72a1)

      return
      end
!------------------------------------------------------------------------------
c
c    text  <  ---  >   integer
c
      subroutine kodir(text,ntext,k,texto,ntexto)
      integer n(80),ntext(80),texto(80)
      character*1 text(80),ntexto(80)
      character*1 nn(80)/'a','b','c','d','e','f','g','h',
     *                 'i','j','k','l','m','n','o','p',
     *                 'q','r','s','t','u','v','w','x',
     *                 'y','z','1','2','3','4','5','6',
     *                 '7','8','9','0','+','|','"','#',
     *                 '.','%',',','''','(',')',' ','=',
     *                 ';','-',':','*','>','<','Z','Y',
     *                 'X','W','V','U','T','S','R','Q',
     *                 'P','O','N','M','L','K','J','I',
     *                 'H','G','F','E','D','C','B','A'/
c     character*1 nn(80)/'a','b','c','d','e','f','g','h',
c    *                 'i','j','k','l','m','n','o','p',
c    *                 'q','r','s','t','u','v','w','x',
c    *                 'y','z','1','2','3','4','5','6',
c    *                 '7','8','9','0','+','|','U','#',
c    *                 'X','Y','Z','''','(',')',' ','=',
c    *                 ';','-',':','*','>','<','?',',',
c    *                 '.','/','W','V','A','B','C','D',
c    *                 'E','F','G','H','I','J','K','L',
c    *                 'M','N','O','P','Q','R','S','T'/
      if(k.eq.1) go to 10
c*** text -> integer    ;    k=0
      do 3 i=1,80
        do 4 j=1,80
           if(text(i).ne.nn(j))go to 5
              ntext(i)=j
              go to 3
    5      continue
    4   continue
    3 continue
      go to 9
c*** integer -> text k=1
   10 continue
      do 13 i=1,80
        do 14 j=1,80
           if(texto(i).ne.j)go to 15
              ntexto(i)=nn(j)
              go to 17
   15      continue
   14   continue
   17 continue
       if(texto(i).lt.1.or.texto(i).gt.80)ntexto(i)=nn(47)
   13 continue
    9 continue
      return
      end
!------------------------------------------------------------------------------
c     . . . Рассчет даты и времени
c
      subroutine dat  (ut,dayt,godt,verno,ut0,ut1,day,god)
      integer god,day,verno,dayt,godt
  900 format(' Ошибка !!! p/p dat:  (ut1-ut0) > utk   !!! '/
     *       ' ut1=',g10.2,'  ut0=',g10.2,'  utk=',g10.2)

      utk=3600*365*24
      if((ut1-ut0).le.utk)go to 4
        print 900,  ut1,ut0,utk
        verno=1
        go to 5
    4 continue
        dayt=day
        godt=god
        ut =ut1
        ng=god
        ngd=ng/4
        ngd=ng-ngd*4
        ngod=365
        if(ngd.eq.0)ngod=366
        i=0
    2   if(ut .lt.86400)go to 1
          ut =ut -86400
          i=i+1
          go to 2
    1   continue
        dayt=dayt+i
        if(dayt.le.ngod)go to 3
          dayt=dayt-ngod
          godt=godt+1
    3   continue
    5 continue
      return
      end
!------------------------------------------------------------------------------
c         Считывание данных Jacobi etc.
      subroutine botread(day,pril,its,ids,kpa,nt) ! add parameter day
                                                  ! day - month day
      dimension pril(kpa,its,ids,nt)
      integer day
      character*2 cday,str*72

      write(cday,'(I2)') day

	print*,'botread cday=',cday
	
      if(day.lt.10)cday(1:1)='0'

      open(9,file='HAMMONIA/PL_'//cday,status='old')
      open(10,file='HAMMONIA/Vwind_'//cday,status='old')
      open(11,file='HAMMONIA/Zwind_'//cday,status='old')
      open(12,file='HAMMONIA/T_'//cday,status='old')
c
     
      do l = 1,nt
       DO i = 1 , ITS
        read(9,*)str

         READ (9,*) (pril(1,i,j,l),j=1,ids)
     !
         read(12,*)str
         READ (12,*) (pril(2,i,j,l),j=1,ids) !Tn
      !
         read(10,*)str
         READ (10,*) (pril(4,i,j,l),j=1,ids) ! meridional
      !
         read(11,*)str
         READ (11,*) (pril(5,i,j,l),j=1,ids)
      enddo
	

      enddo
      do l = 1,nt
       do j=1,ids
        do i=1,its
            pril(3,i,j,l)=0.
            pril(4,i,j,l)=-100.*pril(4,i,j,l)
            pril(5,i,j,l)=100.*pril(5,i,j,l)
        enddo
      enddo

      enddo

  900 format (a25)
  901 format (a1)
  902 format (F5.1,4F8.2)
  903 format (F5.1,3X,E9.2,8X,E9.2,5X,3F8.1)
  904 format (5E11.3)
  101 format(1x,2x,6e12.4)
      
      close(9)
      close(10)
      close(11)
      close(12)
      return
      end