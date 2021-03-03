c    ver     08.11.2017 - var with modeldan
c    ver     21.12.2012 - dtets,ddolgs & dtett,ddolgt - read from danmodel
c            and then compare with inf record  
c            25/05/2012 - pril with kpa & ntIME
c    version 25.03.2012 bmod - modulus of solar magnetic field


      subroutine wwod_mod(mes,god,day,ut0,utk,dtt,dts,tau,
     *          bmpz,bmpy,bmod,vsol,csol,ap0,pkp0,dst0,
     *          ae0,al0,au0,verno,gamma,
     *          ddolgt,dtett,nh,dh,rmin,b,c,ddolgs,dtets,
     *          dlam,rmaxt,nsill,dteta,ns,uprt,uprs,
     *          mas,mast,mass,tet1,tet2,tet3,pdpc,
     *          eps0,om,fac1,fac2,fac3,imja,kpa,ntIME,its,ids)
      
      dimension mas(10),mast(40),mass(30)
    
      integer uprt,uprs,god,day,verno,uth,utm,tkh,tkm
      integer mmes(12)/31,28,31,30,31,30,31,31,30,31,30,31/,b1,c1
      character *8 mes1(12),imja*80
      open(7,file='modeldan',status='old')
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
      if((mas(5).gt.20).or.(mas(4).gt.100)) then
      
         print*,'incorrect mas(4) or mas(5)=', mas(4),mas(5),'*********'
         stop
      end if
      utk=ut0+tk

      read(7,905)bmpy,bmpz,bmod,vsol,csol
      write(10,920)bmpy,bmpz,bmod,vsol,csol
      read(7,907)ap0,pkp0,dst0,ae0,al0,au0
      write(10,924)ap0,pkp0,dst0,ae0,al0,au0


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
  931 format(/2x,f3.0,3x,f3.0,3x,f3.0,2x,d6.2,2x,f4.2,4(2x,d6.2))
      
      write(10,939)tet1,tet2,tet3,eps0,om,fac1,fac2,fac3
  939 format(/'   tet1 = ',f4.0,3x,'tet2 = ',f4.0,5x,'tet3 = ',f4.0/
     *        '   eps0 = ',1pe9.2,7x,'om = ',0pf5.2/
     *        '   fac1 = ',1pe9.2,3x,'fac2 = ',1pe9.2,
     *        '   fac3 = ',1pe9.2)

!  чтение параметров для задания массива PRIL
	read(7,*) 
	

      read(7,*) kpa,ntime

      close(10)
      return
      end


