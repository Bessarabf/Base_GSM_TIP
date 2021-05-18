!  subroutine cyclt_bas, bongi, iacflo_AE, zait, ints, find, inter2, vplm
!     zapvn, copmd, cpf4day, readwr
! function del
c                       CYCLT
C    22/05/2014 ADD PAR2 TO INTERFACE
c    16/04/14 add nl2 to interface
C    25/05/2012 - added kpa & nt - parameters of priL and AE
      subroutine cyclt_bas(god,day,ut0,utk,dtt,dts,tau,b,c,sole,solu,
     *           bmpz,bmpy,bmod,vsol,csol,solen,solet,fa0,fs,ap0,
     *           pkp0,dst0,ae0,al0,au0,nh,ddolgs,dtets,ddolgt,ntsl,
     *           nl,nl2,q,u,dtett,int,ins,rads,ntr,gkoor,kpart,kpars,
     *           nadrt,nadrs,ns,idt,ids,its,nzapt,nzaps,mas,mast,
     *           mass,uprt,kdf,kdu,ldor,isp,rmaxt,par1,PAR2,par,pari,
     *           pole,nr,verno,ni,park,ks,nv,nsu,nse,tet1,tet2,
     *           tet3,pdpc,eps0,om,fac1,fac2,fac3,pglo,
     *           imja,dut,hut,iput,
     *           keut,koob,kut,nara,naou,sut,tut,izap,pril,kpa,nt,
     *           nzapjet,lzapjet,isat)
      character imja*80
      character nara(100)*25,naou(100)*25
      integer god,day,uprt,verno,dayt,godt
      logical readfl
      dimension sole(nse),solu(nsu),solen(nse),solet(nse)
     *          ,mas(10),mast(40),mass(30)
     *          ,ps(10),u(nl),kdf(20),kdu(20),pole(ldor/4),par(nr),
     *          par1(nr),par2(nr),ntsl(nl),q(nv,nl),rads(nh),
     *          dut(isat,100),hut(isat,100),iput(100),keut(100),
     *          kut(100),izap(100),sut(isat,100),tut(isat,100)
     *         ,nzapjet(5),lzapjet(5)
      dimension gkoor(2,its,ids),pari(ni),park(ks),
     *          pril(*)
      dimension pglo(kpars,nh,its,ids)
c
      allocatable vdr(:),potef(:,:,:),
     *          gins(:,:,:,:)
      allocate (vdr(ks),potef(ntr,idt,nl2),
     *          gins(ins,nh,its,ids))
      print *,'cyclt - begin'

	readfl=.true.
      nfile=5
      call zait(readfl,nfile,nzapt,nzaps,nadrt,nadrs,isp,
     *            ldor,kdf,god,day,ut0,ut1,pole)
      uut1=ut1+0.0001
      if(uut1.gt.utk) then
         print 171,ut1,utk
         stop
  171    format('GSMTIP: cyclt - ut1='g10.2,' > utk=',g10.2/
     *     '***** incorrect start and stop time of calculating! ***** ')
      end if
 !     nl2=nl*2+3
      if(mast(23).eq.0) then
c . . . without dreif
        potef= 0.        ! potef(kkk,jjj,iii)=0.
      else
c . . . potential calculate
        readfl=.true.
        n3=ntr*idt*nl2
	  call wpotef(readfl,potef,n3,kdf,kdu,ldor,isp)
      end if
      iqo=mast(25)
      ut=ut0
      nfr=6
      nfw=13
c . . . reading init cond and writing in iteratioon part of tube
      call copmd(nfr,nfw,kdf,isp,ldor,kdu,par,nr,mass,mast)

      do i=1,ks
        vdr(i)=0.
      end do
      nomsl=1
      md=1
      dolm=0.
      readfl=.true.
      kpar=kpart
      nfile=6
c . . . reading init cond and writing interpolation parameters
      do i=1,idt
        call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
     *           kdf,ldor,isp,md,par,pole,nr,mast)
        call ints(dolm,par,nr,rads,nh,ni,pari,ins,its,park,ks,
     *            ntsl,nl,ntr,kdf,ldor,isp,pole,kpart,vdr,dtett,
     *            u,ddolgs,dtets,gins,ids)
        dolm=dolm+ddolgt
      end do
c . . . data in poles
      call bongi(gins,ins,nh,its,ids)
c     . . . continue or new calculation?
      if(uprt.eq.3) then
        call dat(ut,dayt,godt,verno,ut0,ut1,day,god)
      else
        ut1=ut0
        dayt=day
        godt=god
c . . . date and time write in inf record
        read(4,rec=1) pole
        pole(2)=god
        pole(3)=day
        pole(4)=ut0
        pole(5)=ut1
        write(4,rec=1) pole
c . . .
      end if
      delta=del(dayt)
      utn=ut1+tau
      md=0
      dtt0=dtt
      if(dtt0.gt.dts) dtt0=dts
!!!!!!!!!!!! cycle on time begin  !!!!!!!!!!!!!!!!!!!!! 
  121 if(utn.ge.(utk+tau))go to 21
        uts=ut1+dts
        if(uts.ge.86400.) then
          call dat(ut,dayt,godt,verno,ut0,uts,day,god)
          delta=del(dayt)
        end if
  122   if(uts.ge.utn)go to 22
!!!!!!!!!!!! from danmodel
          pkpj1=pkp0
          AEpdpc=ae0
          AEj2=ae0
!!!!!!!!!!!!!!!!!!!!!!!!!!
          call iacflo_AE(uts,sole,solu,bmpz,bmpy,mas,
     *               csol,vsol,fa0,pkp0,ap0,ae0,dst0,al0,au0,
     *               fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps,
     *               AEj2)
          call cyclt1_bas(god,day,ut0,ut1,dtt,dts,tau,solen,sole,solu,
     *              solet,bmpz,bmpy,bmod,vsol,csol,ntr,gkoor,ps,gins,
     *              fa,fs,ap,pkp,dst,ae,al,au,rads,ni,nv,vdr,nh,
     *              ddolgs,dtets,ddolgt,ntsl,nl,idt,ids,its,q,u,
     *              kpart,kpars,nadrt,nadrs,nzapt,nzaps,mast,mass,
     *              mas,dtt0,uts,ns,dtett,int,ins,b,c,kdf,kdu,ldor,
     *             isp,rmaxt,nsu,nse,par1,PAR2,par,park, ks ,pari,pole,
     *              nr,verno,delta,dayt,godt,tet1,tet2,tet3,pdpc,
     *              eps0,om,fac1,fac2,fac3,nl2,potef,iqo,pglo,imja,
     *              dut,hut,iput,keut,koob,kut,nara,naou,sut,tut,izap,
     *              pril,kpa,nt,nzapjet,lzapjet,AEpdpc,AEj2,pkpj1,isat)

          if(verno.eq.1)uts=utn
          uts=uts+dts
          if(uts.ge.86400.) then ! ut hour > 24 h, new day
            call dat(ut,dayt,godt,verno,ut0,uts,day,god)
            delta=del(dayt)
          end if
          if(uts.ge.utn) then
            dts=utn-uts+dts
            uts=utn
            if(uts.ge.86400.) then
              call dat(ut,dayt,godt,verno,ut0,uts,day,god)
              delta=del(dayt)
            end if
          end if
          go to 122
   22   continue
!!!!!!!!!!!! from danmodel
          pkpj1=pkp0
          AEpdpc=ae0
          AEj2=ae0
!!!!!!!!!!!!!!!!!!!!!!!!!!
        call iacflo_AE(uts,sole,solu,bmpz,bmpy,mas,
     *             csol,vsol,fa0,pkp0,ap0,ae0,dst0,al0,au0,
     *               fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps,
     *               AEj2)
        call cyclt1_bas(god,day,ut0,ut1,dtt,dts,tau,solen,sole,solu,
     *              solet,bmpz,bmpy,bmod,vsol,csol,ntr,gkoor,ps,gins,
     *              fa,fs,ap,pkp,dst,ae,al,au,rads,ni,nv,vdr,nh,
     *              ddolgs,dtets,ddolgt,ntsl,nl,idt,ids,its,q,u,
     *              kpart,kpars,nadrt,nadrs,nzapt,nzaps,mast,mass,
     *              mas,dtt0,uts,ns,dtett,int,ins,b,c,kdf,kdu,ldor,
     *             isp,rmaxt,nsu,nse,par1,PAR2,par,park,ks,pari,pole,
     *              nr,verno,delta,dayt,godt,tet1,tet2,tet3,pdpc,
     *              eps0,om,fac1,fac2,fac3,nl2,potef,iqo,pglo,imja,
     *              dut,hut,iput,keut,koob,kut,nara,naou,sut,tut,izap,
     *              pril,kpa,nt,nzapjet,lzapjet,AEpdpc,AEj2,pkpj1,isat)
        ut1=utn
        utn=utn+tau
        readfl=.false.
        nfile=7
        nrazm=ins*nh*its*ids
        nkd=20
        nf=4
c . . . end time-step and writing into f4
        call readwr(readfl,nfile,kdf,kdu,nkd,nrazm,gins,ldor,isp,nf)
c . . . copy f4 or not?
        if(mas(9).ne.0) then
          iuts=uts+0.0001
          idts=mas(9)*3600+0.0001
          if(mod(iuts,idts).eq.0) then
            isf4=0
            do i=1,8
               isf4=isf4+kdu(i)
            end do
!            call copyf4(uts)
!            call copyf4day(uts,day)
             call cpf4day(uts,day,isf4,ldor)
          end if
        end if
	go to 121
!!!!!!!!!!!!!!!!! end of time step  !!!!!!!!!!!!!!!!!!!!!!!!

   21 continue   !!!!!!!!!!!!!!!!!!!  end of calculating time 
c . . . writing time at the end of program
        call dat(ut,dayt,godt,verno,ut0,ut1,day,god)
        readfl=.false.
        call zait(readfl,5,nzapt,nzaps,nadrt,nadrs,isp,
     *           ldor,kdf,god,day,ut0,ut1,pole)
        call zait(readfl,6,nzapt,nzaps,nadrt,nadrs,isp,
     *           ldor,kdf,god,day,ut0,ut1,pole)
c . . . writing time and date into inf record
      read(4,rec=1) pole
      pole(2)=god
      pole(3)=day
      pole(4)=ut0
      pole(5)=ut1
      write(4,rec=1) pole
c . . .
c  26 continue
      print *,'GSMTIP: cyclt - end'
      deallocate (vdr,potef,gins)
      return
      end
!------------------------------------------------------------------------------
      subroutine bongi(gins,ins,nh,its,ids)
      dimension gins(ins,nh,its,ids)
      i2=its-1
      do 3 i=1,3
          np=i
          if(i.eq.2)np=5
          if(i.eq.3)np=6
       do 2 k = 1 , nh
c      ssp - sum s.pole
c      snp - sum n.pole
        s np=0.
        s sp=0.
        do 1 j = 1 , ids
         snp=snp+gins(np,k,2,j)
         ssp=ssp+gins(np,k,i2,j)
   1    continue
c
        unp=snp/ids
        usp=ssp/ids
        do 4 j=1,ids
          gins(np,k,1,j)=unp
          gins(np,k,its,j)=usp
    4   continue
    2  continue
    3 continue
      return
      end
!------------------------------------------------------------------------------
      function del(day)
      integer day
      data tg/4.3481234e-1/,p/1.7214206e-2/
      del=atan(tg*sin(p*(day-80.)))
      return
      end
!------------------------------------------------------------------------------
c     . . . energy and flux of electron precipitations
      subroutine iacflo_AE(utn,sole,solu,bmpz,bmpy,mas,
     *                  csol,vsol,fa0,pkp0,ap0,ae0,dst0,al0,au0,
     *                  fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps,
     *                  AEj2)
      dimension sole(nse),solu(nsu),mas(10),solen(nse),ps(10)
! . . . soft electron (cusp) 70 deg
!     ps(1 )=1.e8
      ps(1 )=5.e8
!     ps(1 )=5.e9
!     ps(1 )=2.5e9
!     ps(1 )=1.5e9
      ps(2 )=1.
      ps(3 )=0.20e3
! . . . hard electron (auroral) 70 deg
!     ps(4 )=5.e8
      ps(4 )=1.e7  ! flux
      ps(5 )=1.0   ! param gamma
      ps(6 )=3.e3  ! characteristic energy
!     ps(6 )=5.e3
! . . .soft electrons 80 deg
      ps(7 )=0.e8  ! first zone flux
!     ps(7 )=5.e8
!     ps(7 )=4.e9
      ps(8 )=0.1e3 ! characteristic energy
!     ps(8 )=0.05e3
      ps(9 )=0.e8  ! second zone flux
!     ps(9 )=3.e8
      ps(10)=0.05e3 ! characteristic energy
! . . .geomagnetic indexes ...
      if(mas(2).eq.0) then
!       . . .from input data
        fa=fa0
        pkp=pkp0
        ap=ap0
        ae=ae0
        dst=dst0
        al=al0
        au=au0
!       . . .calculating
      else
        fa=1
        pkp=1
        ap=1
        ae=1
        dst=1
        al=1
        au=1
      end if
!!!!!!!!!!!     storm ..............
!    . . . «ј¬»—»ћќ—“№ ¬џ—џѕјЌ»… ќ“ Kp
!1
!         ps(4)=ps(4)*(-3.5+6.5*pkpj2)
!1
!2
!         ps(4)=ps(4)*1.43*pkpj2
!2
!3
         ps(4)=ps(4)*(1.0+0.004*AE0)
!3
!         ps(6)=5.e3
         ps(6)=2.58e3+0.003e3*AE0
!!!!!!!!!!!     storm ..............
      return
      end
!------------------------------------------------------------------------------
      subroutine wpotef(readfl,potef,n,kdf,kdu,
     *                   ldor,isp)
      dimension potef(n)
	allocatable pole(:)

      dimension kdf(20),kdu(20)
      logical readfl

	allocate (pole(ldor/4))

      mdor=ldor/4
      nfile=4
      nf=5
      if(nfile.le.8)nf=4
      isp=kdf(nfile)+1
      jj=kdu(4)
      jm=jj-1
      l=jj*mdor
      lm=l-mdor
      ll=n-lm
      if(readfl) go to 4
        k=1
        do 2 i=1,jm
          do 1 j=1,mdor
            pole(j)=potef(k)
            k=k+1
    1     continue
          write(nf,rec=isp)pole
          isp=isp+1
    2   continue
        do 3 j=1,ll
          pole(j)=potef(k)
          k=k+1
    3   continue
        write(nf,rec=isp)pole
        go to 8
    4 continue
        k=1
        do 6 i=1,jm
          read(nf,rec=isp)pole
          do 5 j=1,mdor
            potef(k)=pole(j)
            k=k+1
    5     continue
          isp=isp+1
    6   continue
        read(nf,rec=isp)pole
        do 7 j=1,ll
          potef(k)=pole(j)
          k=k+1
    7   continue
    8 continue
      deallocate (pole)
      return
      end
!------------------------------------------------------------------------------
c . . . «апись и чтение времени и даты
      subroutine zait(readfl,nfile,nzapt,nzaps,nadrt,nadrs,isp,
     *                ldor,kdf,god,day,ut0,ut1,pole)
      integer nfa(20)/4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,
     *       0,0,0,0,0/
      dimension kdf(20),pole(ldor/4)
      logical readfl
      integer adres,god,day
      mdor=ldor/4
      nzap=nzapt
      adres=nadrt
      if(nfile.eq.6)go to 4
        nzap=nzaps
        adres=nadrs
    4 continue
      isp=kdf(nfile)+nzap
      nf=nfa(nfile)
      read(nf,rec=isp)(pole(i),i=1,mdor)
      if(readfl) go to 1
        pole(adres+3)=god
        pole(adres+2)=day
        pole(adres+1)=ut1
        pole(adres)=ut0
        isp=kdf(nfile)+nzap
        write(nf,rec=isp)(pole(i),i=1,mdor)
        go to 2
    1 continue
        ut1=pole(adres+1)
    2 continue
c      print 900,ut0,ut1,god,day,readfl
      return
      end
!------------------------------------------------------------------------------
      subroutine  ints(dolm,par,nr,rads,nh,ni,pari,ins,its,park,
     *                 ks,ntsl,nl,ntr,kdf,ldor,isp,pole,kpart,vdr,
     *                 dtett,u,ddolgs,dtets,gins,ids)

      dimension ntsl(nl),par(nr),rads(nh),park(ks),gins(ins,nh,its,ids),
     *      pari(ins,nh,its),kdf(20),pole(ldor/4),vdr(ks),
     *      msum(45),u(nl)
      allocatable vn1(:),vn2(:),plm(:),grz(:,:)
      double precision pi,f,re,cr
      logical readfl
      data  pi/3.14159265359d0/
	allocate (vn1(ins),vn2(ins),plm(ins),grz(ins,nh))
  900 format(' ',6e12.4)
      cr=180.d0/pi
	!!!
	pari=0.
	!!!
      
      msum(1)=0
      do  nomsl=2,nl
        msum(nomsl)=msum(nomsl-1)+ntsl(nomsl-1)
      end do
      re=6.37102d8
C  Ёкватор: нижний узел вычисл€етс€ интерпол€цией
C            из NL-й линии по горизонтали
      i=msum(nl)
      call zapvn(i,kpart,par,vdr,nr,ks,vn1,ins)
      x1=park(i*2+2)
      i=i+ntsl(nl)-1
      call zapvn(i,kpart,par,vdr,nr,ks,vn2,ins)
      x2=park(i*2+2)
      xi=90.
C  »нтерпол€ци€:
      call inter2(x1,x2,xi,vn1,vn2,plm)
      do  i=1,ins
        grz(i,ntr)=plm(i)
      end do
      tet=90.
      p4=plm(3)
      p3=plm(4)
C  ѕересчет составл€ющих скорости:
      call vplm(plm(2),p3,tet)
      plm(3)=p3
      plm(4)=p4
      k=its/2+1
      do i=1,ins
        pari(i,ntr,k)=plm(i)
      end do
C  ќстальные точки экватора:
      ro=re/(re+rads(ntr))
      tet1=180.-park(2)
      nt1=ntr+1
      do 7 nt=nt1,nh
        unt=re/(re+rads(nt))
C       TT - коширота точки пересечени€ с NTR-й высотой силовой линии
C       с вершиной на NT-й высоте, град.
        tt=sqrt(unt/ro)
        tt=asin(tt)
        tt=tt*cr
        n=(tt-tet1)/dtett+1
        i=msum(n)+ntsl(n)/2
        call zapvn(i,kpart,par,vdr,nr,ks,vn2,ins)
        x2=park(i*2+1)
        xi=rads(nt)
        if(n.ne.nl) go to 8
          do i=1,ins
            vn1(i)=grz(i,ntr)
          end do
          x1=rads(ntr)
          go to 10
    8   continue
          i=msum(n+1)+ntsl(n+1)/2
          call zapvn(i,kpart,par,vdr,nr,ks,vn1,ins)
          x1=park(i*2+1)
   10   continue
        call inter2(x1,x2,xi,vn1,vn2,plm)
        do i=1,ins
          grz(i,nt)=plm(i)
        end do
        p4=plm(3)
        p3=plm(4)
c       call vplm(plm(2),plm(3),tet)
ccc     plm(2)=0.
        call vplm(plm(2),p3,tet)
        plm(3)=p3
        plm(4)=p4
        do i=1,ins
          pari(i,nt,k)=plm(i)
        end do
    7 continue
C   онец расчета на экваторе
C Ќачало цикла по коширотам от северного полюса к южному:
      k=2
      tet=dtets
   20 if(tet.ge.180.d0) go to 13
        do nt=ntr,nh
          h=rads(nt)
          t=tet/cr
          f=sin(t)
          !unt=re/(re+h)*(sin(t))**2
          unt=re/(re+h)*f*f
          
          call find(nl,unt,u,n)
          if(n.eq.nl)then
            ntt=nt-ntr
            if(tet.lt.90.)ntt=ntsl(n)-ntt-1
            i=msum(n)+ntt
            i1=i*kpart+1
            i2=i*2+1
            plm(1)=par(i1)
            plm(2)=par(i1+3)
            plm(3)=vdr(i2)
            plm(4)=vdr(i2+1)
            plm(5)=par(i1+7)
            plm(6)=par(i1+6)
          else
            hsl=re/u(n+1)-re
            ntt=nt-ntr
            if(tet.lt.90.)ntt=ntsl(n)-ntt-1
            i=msum(n)+ntt
            call zapvn(i,kpart,par,vdr,nr,ks,vn1,ins)
            x1=park(i*2+2)
            if(h.le.hsl) go to 15
              do i=1,ins
                vn2(i)=grz(i,nt)
              end do
              x2=90.
              go to 17
   15       continue
              ntt=nt-ntr
              if(tet.lt.90.)ntt=ntsl(n+1)-ntt-1
              i=msum(n+1)+ntt
              call zapvn(i,kpart,par,vdr,nr,ks,vn2,ins)
              x2=park(i*2+2)
   17       continue
            xi=tet
            call inter2(x1,x2,xi,vn1,vn2,plm)
          end if
          p4=plm(3)
          p3=plm(4)
c         call vplm(plm(2),plm(3),tet)
ccc       plm(2)=0.
          call vplm(plm(2),p3,tet)
          plm(3)=p3
          plm(4)=p4
          do i=1,ins
            pari(i,nt,k)=plm(i)
          end do
        end do
        tet=tet+dtets
        k=k+1
        if(abs(tet-90.).gt.0.001)go to 19
          tet=tet+dtets
          k=k+1
   19   continue
        go to 20
   13 continue

      m=dolm/ddolgs+1
      do 23 i=1,ins
        do 22 j=1,nh
          do 21 k=1,its
            !if(ISNAN(pari(i,j,k))) then
            !  print*, 'ins=',i,' alt=',rads(j),' its=',k,pari(i,j,k),
     *      ! pari(i,j-1,k), 'dolg=',m
            !  pause
            !end if
            gins(i,j,k,m)=pari(i,j,k)
   21     continue
   22   continue
   23 continue
      deallocate (vn1,vn2,plm,grz)
      return
      end
!------------------------------------------------------------------------------
      subroutine find(n,u,x,m)
      dimension x(n)
      data i/1/
c     if(u.gt.x(n).or.u.le. x(1)) print 900,u,x(1),x(n)
c 900 format(' find'/
c    *' значение u выходит за пределы массива x в ',a4/
c    *'   u=',e10.3,'  x(1)=',e10.3,'  x(n)=',e10.3)
      if(i.ge.n) i=1
      if(u.lt.x(i)) go to 10
      if(u.le.x(i+1)) go to 30
c
   10 i=1
      j=n+1
   20 k=(i+j)/2
      if(u.lt.x(k)) j=k
      if(u.ge.x(k)) i=k
      if(j.gt.i+1) go to 20
   30 m=i
      return
      end
!------------------------------------------------------------------------------
      subroutine inter2(x1,x2,xi,vn1,vn2,plm)
      dimension plm(6),vn1(6),vn2(6)
      a1=amin1(x1,x2)
      a2=amax1(x1,x2)
      ! 10.03.2018 first and lust tube is vertical
      if(xi.lt.a1.and.a1.eq.x1.or.xi.gt.a2.and.a2.eq.x1)then
        do i=1,6
          plm(i)=vn1(i)
        end do
      end if
      if(xi.lt.a1.and.a1.eq.x2.or.xi.gt.a2.and.a2.eq.x2)then
        do i=1,6
          plm(i)=vn2(i)
        end do
      end if
      if(xi.ge.a1.and.xi.le.a2)then
        hi=x2-x1
        w1=(xi-x1)/hi
        if(xi.lt.a1.or.xi.gt.a2)print*,' ints x1,xi,x2  ',x1,xi,x2
        w2=1.-w1
c       if(vn1(1).le.0.)vn1(1)=1.e15
c       if(vn2(1).le.0.)vn2(1)=1.e15
c       p1=alog(vn1(1))
c       p2=alog(vn2(1))
c       plm(1)=exp(w1*p2+w2*p1)
        do i=1,6
          plm(i)=w1*vn2(i)+w2*vn1(i)
        end do
      end if
      return
      end
!------------------------------------------------------------------------------
      subroutine vplm(p1,p2,tet)
      double precision f,tr,sn,cs,tn2,si
      f=0.01745329252
      tr=(90.-tet)*f
      sn=dsin(tr)
      cs=dcos(tr)
      tn2=(sn/cs)*2
      si=datan(tn2)
      sn=dsin(si)
      cs=dcos(si)
      a1=p1
      a2=p2
      p1=-(a1*sn+a2*cs)
      p2=a2*sn-a1*cs
      return
      end
!------------------------------------------------------------------------------
c
      subroutine zapvn(i,kpart,par,vdr,nr,ks,vn,ins)
      dimension par(nr),vdr(ks),vn(ins)
      i1=i*kpart+1
      i2=i*2+1
      vn(1)=par(i1)
      vn(2)=par(i1+3)
      vn(3)=vdr(i2)
      vn(4)=vdr(i2+1)
      vn(5)=par(i1+7)
      vn(6)=par(i1+6)
      return
      end
!------------------------------------------------------------------------------
      subroutine copmd(nfr,nfw,kdf,isp,ldor,kdu,par,nr,mass,mast)
      dimension par(nr),mast(40),mass(30),kdu(20),kdf(20)
  900 format(' copmd:',20i5)
  901 format(' copmd:  ********  error *********'/
     *       ' размеры участков на дисках',
     *       ' не совпадают'/
     *  ' kdu(',i2,')=',i7,'   kdu(',i2,')=',i7,' !!!!!!! STOP !')
  902 format(' copmd:  ********  error ********** '/
     *       ' размерность области I/O меньше mdor'/
     *         ' nr = ',i10,'  <   mdor=',i5,'  !!!!!! STOP  !')
  903 format(' p/p copmd: nfr=',i5,'  ---> ','  nfw=',i5)
      
	
	mdor=ldor/4
      if(nfr.eq.14.and.mast(13).eq.0)go to 9
        if(nfw.eq.14.and.mast(13).eq.0)nfw=13
        
        print 900,nfr,nfw
        if(nr.lt.mdor) go to 3
          i1=kdf(nfr)+1
      	i2=kdf(nfw)+1
          kd1=kdu(nfr)
          kd2=kdu(nfw)
          if(kd1.ne.kd2) go to 2
            nf1=5
            nf2=5
            if(nfr.le.8) nf1=4
            if(nfw.le.8) nf2=4
            do 1 j=1,kd1
              read(nf1,rec=i1)(par (i),i=1,mdor)
              write(nf2,rec=i2)(par (i),i=1,mdor)
              i1=i1+1
              i2=i2+1
    1       continue
            print 903,nfr,nfw
            go to 5
    2     continue
            print 901,nfr,kd1,nfw,kd2
            stop
    5     continue
          go to 8
    3   continue
          print 902,nr,mdor
          stop
    8   continue
        print *,' copmd - end'
    9 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine cpf4day(uts,day,nzap,ldor)
      allocatable pole(:)
      integer day
      character sr*2,srd*3
      character*80 fname
      allocate( pole(ldor/4))
      ihour=uts/3600.+0.001
	write(sr,'(i2)') ihour
        write(srd,'(i3)') day
	!print*,sr
	if(ihour.lt.10) sr(1:1)='0'
      if(day.lt.100) srd(1:1)='0'
      if(day.lt.10) srd(1:2)='00'
      fname='file4_'//srd//'.'//sr  
      open(44,file=(fname),
     *       access='direct',form='unformatted',recl=ldor) 
        do irecl=1,nzap
   !        print*,'irecl=',irecl
           read(4,rec=irecl) pole
           write(44,rec=irecl) pole
        end do
        close(44)
        deallocate(pole)
	return
	end
!------------------------------------------------------------------------------
c
c     READWR.f     11.03.93 Naumowa
c    cyclt:readwr
c    cyclt:cyclt1:readwr
c    cyclt:cyclt1:sklad:readwr
      subroutine readwr (readfl,nfile,kdf,kdu,nkd,nrazm,param,
     *                  ldor,isp,nf)
      dimension param(nrazm),kdu(nkd),kdf(nkd)
      logical readfl
      print *,' readwr: gins --> disk'
      print *,' readwr:  nf=',nf,'   nfile=',nfile
      mdor=ldor/4
      nk=mdor*kdu(nfile)
      if(nk.lt.nrazm)go to 1
        isp=kdf(nfile)+1
c       nf=4
c       if(nfile.gt.8)nf=5
        ind=1
        indk=mdor
        if(nrazm.lt.mdor)go to 8
          ncycl=nrazm/mdor
          nost=nrazm-ncycl*mdor
          if(nost.gt.0)ncycl=ncycl+1
          go to 9
    8   continue
          indk=nrazm
          ncycl=1
    9   continue
        if(readfl)go to 2
          do 5 i=1,ncycl
            write(nf,rec=isp)(param(j),j=ind,indk)
            isp=isp+1
            ind=indk+1
            indk=indk+mdor
            if(indk.gt.nrazm)indk=nrazm
    5     continue
          go to 3
    2   continue
          do 6 i=1,ncycl
            read(nf,rec=isp)(param(j),j=ind,indk)
            isp=isp+1
            ind=indk+1
            indk=indk+mdor
            if(indk.gt.nrazm)indk=nrazm
    6     continue
    3   continue
        go to 4
    1 continue
        print 900,nrazm,nfile,kdu(nfile),ldor
  900 format(' p/p readwr:  ****** ERROR !! *****'/
     *       '   nrazm=',i5,'      kdu(',i3,') =',i10)
        stop
    4 continue
      return
      end
!------------------------------------------------------------------------------

