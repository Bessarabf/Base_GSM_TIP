c                       CYCLT
c    23/03/2019 Ionize with HAMMONIA ionization 
c    18.05.18 Ion Drag sent to HAMMONIA
c    10/05/2018 Joul heating sent to HAMMONIA
C    22/05/2014 ADD PAR2 TO INTERFACE
c    16/04/14 add nl2 to interface
C    25/05/2012 - added kpa & nt - parameters of priL and AE
      subroutine cyclt_HAM(god,day,ut0,utk,dtt,dts,tau,b,c,sole,solu,
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
          call cyclt1_HAM(god,day,ut0,ut1,dtt,dts,tau,solen,sole,solu,
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
        call cyclt1_HAM(god,day,ut0,ut1,dtt,dts,tau,solen,sole,solu,
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
