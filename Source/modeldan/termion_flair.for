! Новая версия модели, адаптированная под любую сетку
!
! НИКАКИХ INCLUDE - файлов!
!
! Параметры файла f4 записаны в информационной записи
! отсюда же используем информацию о сетке и шагах 
!
! Максимально используем динамическое распределние массивов

! объявление статических массивов
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Parameter(isat=8700,ldor0=4096)
!	! Parameter(kpa=5,ntIME=4)! number of the parameters and time moments for PRIL Liu_dat
!	Parameter(kpa=5,ntIME=12) ! number of the parameters and time moments for PRIL Jacobi_dat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character nara(100)*25,naou(100)*25
      character imja*80

      integer uprt,uprs,verno/0/,god,day
	dimension mas(10),mast(40),mass(30),kdf(20),kdu(20),pole(ldor0/4)

      dimension
     *          
     *          dut(isat,100),hut(isat,100),iput(100),
     *          keut(100),kut(100), 

     *          sut(isat,100),tut(isat,100),izap(100),
     *          nzapjet(5),lzapjet(5)
! объявление динамических массивов
      allocatable
     *            q(:),gkoor(:)
     *           ,park(:),pari(:),par(:),par1(:),par2(:)
     *           ,pglo(:)
     *           ,pril(:,:,:,:) ! for other lowboundary model different numbers
     *           ,sole(:),solu(:),solen(:),solet(:)
     *           ,ntsl(:),u(:),rads(:)
 
	data  lzapjet/20,20,12,12,24/,nzapjet/5*0/

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

!
	
! чтение modeldan
      call wwod_mod(mes,god,day,ut0,utk,dtt,dts,tau,
     *          bmpz,bmpy,bmod,vsol,csol,ap0,pkp0,dst0,ae0,
     *          al0,au0,verno,gamma,
     *          ddolgt,dtett,nh,dh,rmin,b,c,ddolgs,dtets,
     *          dlam,rmaxt,nsill,dteta,ns,uprt,uprs,
     *          mas,mast,mass,tet1,tet2,tet3,pdpc,
     *          eps0,om,fac1,fac2,fac3,imja,kpa,ntIME,its,ids)
!
      nse=mas(4)
      nsu=2*mas(5)
      allocate(sole(nse),solu(nsu),solen(nse),solet(nse))
      call solusolen(solu,sole,solen,solet,nsu,nse,mas)     
      print*,sole
      print*,solu
      print *, 'MES=', mes

! сравнение на непротиворечие входных данных из информационной записи и danmodel
      if(verno.NE.0) then
              print*,' WWOD__D: incorrect input data between',
     *        'danmodel and INF record! ERRCODE=',verno
	          stop
      end if

! выделение памяти под динамические массивы
      ni=ntpl*int ! размер массива трубки (плоскость)
      ks=ntpl*2
      nr=nh*its*kpars	! размер массива шара (плоскость)
   
      allocate(pglo(kpars*nh*its*ids),
     *         par(nr),par1(nr),par2(nr),park(ks),pari(ni),
     *         q(nv*nl),gkoor(2*its*ids),
     *         pril(kpa,its,ids,ntIME),
     *	       rads(nh),u(nl) )

! расчетная часть
!  чтение файла fiza
	  if(mast(31).ne.0) then 
         call rpco(dut,hut,iput,keut,koob,kut,nara,naou,sut,tut,
     *             isat)

         open(9,file='fiza',status='old')
         read(9,'(10i5)')izap
         close(9)
	  end if
  !    nse=mas(4)
  !    nsu=mas(5)*2


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
        !call botread(mes,day,pril,its,ids,kpa,ntIME)  !For Hagan
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
	
! освобождение памяти
      deallocate (pglo,par,par1,PAR2,q,gkoor,pril,
     *            sole,solu,solen,solet,rads,u,ntsl)  

! STOP
      stop
      end