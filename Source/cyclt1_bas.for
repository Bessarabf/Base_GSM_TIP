cc                       CYCLT1
c   _bas - base variant based on EAGLE 

C            22/05/2014 - ADD PAR2 TO INTERFACE
c            25/05/2012 - kpa & nt to priL
c    version 25.03.2012 - bmod and isat
c    
      subroutine cyclt1_bas(god,day,ut0,ut1,dtt,dts,tau,solen,sole,solu,
     *                 solet,bmpz,bmpy,bmod,vsol,csol,ntr,gkoor,ps,
     *                 gins,fa,fs,ap,pkp0,dst,ae,al,au,rads,ni,nv,vdr,
     *                 nh,ddolgs,dtets,ddolgt,ntsl,nl,idt,ids,its,q,u,
     *                 kpart,kpars,nadrt,nadrs,nzapt,nzaps,mast,mass,
     *                 mas,dtt0,uts,ns,dtett,int,ins,b,c,kdf,kdu,ldor,
     *                isp,rmaxt,nsu,nse,par1,PAR2,par,park,ks,pari,pole,
     *                 nr,verno,delta,dayt,godt,tet1,tet2,tet3,pdpc,
     *                 eps0,om,fac1,fac2,fac3,nl2,potef,iqo,pglo,
     *                 imja,dut,hut,iput,keut,
     *                 koob,kut,nara,naou,sut,tut,izap,pril,kpa,nt,
     *                 nzapjet,lzapjet,AEpdpc,AEj2,pkp,isat)
!
   
      character *1 imja(80)
      character *25 namefl
      dimension par(nr),gkoor(2,its,ids),ps(10),
     *          pari(ni),par1(nr),par2(nr),vdr(ks),park(ks),
     *          potef(ntr,idt,nl2),pglo(kpars,nh,its,ids),
     *          gins(ins,nh,its,ids),
     *          sole(nse),solu(nsu),solet(nse),pole(ldor/4),solen(nse),
     *          nzapjet(5),lzapjet(5),
     *          ntsl(nl),q(nv,nl),u(nl),
     *          mast(40),mass(30),kdf(20),kdu(20),rads(nh),mas(10),
     *          dut(isat,100),hut(isat,100),iput(100),keut(100),
     *          kut(100),sut(isat,100),tut(isat,100),izap(100),pril(*)
      
      allocatable E0(:,:),FAE(:,:),vert(:,:),qom(:)
      

      character nara(100)*25,naou(100)*25
      integer god,day,verno,godt,dayt
      logical readfl
      !data pi/3.14159265359/
      print *,'GSMTIP: cyclt1 - begin'
      CALL CPU_TIME(TT1)
      allocate (E0(its,ids),FAE(its,ids),vert(kpart,nv),qom(nl))
      i=0
!!!!  do j=1,nl
        qom=0.
!!!!  end do
      qmax=0.
!!!!  pdpc cross potential relations
!!!!  pdpc(Kp) with Zhang&Pacston precipitations mass(6)=2 or 6
!!!!  pdpc(Ae) with Vorobiev and Yagodkina 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mass(12).eq.2.or.mass(12).eq.6) then
        pdpc=26.4+13.3*pkp0
      else if(mass(12).eq.3) then
        pdpc=38.0+0.089*AEpdpc
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    shift tet2 on latitude with dependences with pdpc 
      if(pdpc.ge.0..and.pdpc.lt.40.)tet2=25.
      if(pdpc.ge.40..and.pdpc.lt.50.)tet2=30.
      if(pdpc.ge.50..and.pdpc.lt.88.5)tet2=35.
      if(pdpc.ge.88.5.and.pdpc.lt.127.)tet2=40.
      if(pdpc.ge.127..and.pdpc.lt.165.4)tet2=45.
      if(pdpc.ge.165.4.and.pdpc.lt.200.)tet2=50.
      if(pdpc.ge.200.)tet2=55.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      fac2=3.06e-13+2.14e-16*ae
!      fac2=(3.e-14+0.006*1.e-14*ae)*1.28
!      fac2=(3.e-14+0.006*1.e-14*AEpdpc)*1.28
!      fac2=(3.e-13+0.006*1.e-13*AEpdpc)*1.28
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       fac2 from danmodel only !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    1 if(i.ge.ns) go to 2
        dtt=dtt0
        utt=uts+dtt-dts
        alt=rads(ntr)
!!!   electron precipitation
      ! Zhang & Paxton
      if(mass(12).eq.2.or.mass(12).eq.6)then
        call aepp(mass,pkp,nl2,idt,ddolgs,dtets,uts,delta,
     *            E0,FAE,its,ids)
      ! Vorobiev & Yagodkina
      else if(mass(12).eq.3) then
        call aeppVY(mass,nl2,idt,ddolgs,dtets,al,dst,uts,delta,E0,FAE)
      else  
        print*, 'incorrect mass(12) - var of precipitations'
      end if
!------------------------------------------------------------------------------
      call shar_bas(day,god,dayt,godt,uts,tau,dts,sole,solu,solen,nsu,
     *           nse,kpars,ins,int,rads,nh,gkoor,its,ids,ddolgs,
     *           dtets,dtett,fa,fs,ap,pkp,dst,ae,al,au,bmpz,bmpy,
     *           vsol,ps,csol,mass,delta,kdu,kdf,ldor,isp,pole,par,
     *           nr,pari,par1,ni,ddolgt,ntsl,nl,verno,park,ks,potef,
     *           nl2,ntr,gins,solet,ut0,qom,qmax,iqo,mast,pglo,
     *           pril,kpa,nt,E0,FAE)
!------------------------------------------------------------------------------
      do ix = 1 , kpars
        print *,'pglo (',ix,') max=',maxval(pglo(ix,:,:,:)),'min=',
     *    minval(pglo(ix,:,:,:))
      end do
!------------------------------------------------------------------------------
!     Potencial calculation 
       if(mast(13).eq.0) go to 8
          if(mast(18).eq.1) go to 119
                  potef=0.
  119     continue
           if(mast(26).eq.0)then
            print *,' cyclt1:  potential - Volland'
            call potvol(delta,uts,pkp,bmpy,alt,ddolgt,dtett,
     *           potef,idt,nl2,ntr,b,c,isp,readfl,kdf,kdu,ldor)
          else if(mast(26).eq.1)then
              print*,' cyclt1: potential integrate on height'
cc
              dekr=1.
              call potnew(delta,uts,pkp,bmpy,bmpz,vsol,csol,ntr,ap,dst,
     *                al,ae,au,rads,nh,ddolgs,ddolgt,dtets,dtett,
     *                isp,readfl,kpars,kdf,kdu,nr,ldor,nfile,md,idt,
     *                its,nl2,par,pole,potef,mast,tet1,tet2,tet3,
     *                pdpc,eps0,om,fac1*dekr,fac2*dekr,fac3*dekr,
     *                pglo,ids,nzapjet,lzapjet)
          else if(mast(26).eq.2)then
                print *,' cyclt1: potential integrate along',
     *			  ' a geomagnetic field - 1st var'
                call potnewmf1(delta,uts,pkp,bmpy,bmpz,
     *			 vsol,csol,ntr,ap,
     *             dst,al,ae,au,rads,nh,ddolgs,ddolgt,dtets,dtett,
     *             isp,readfl,kpars,kdf,kdu,nr,ldor,nfile,md,idt,
     *             its,nl2,par,pole,potef,mast,tet1,tet2,tet3,
     *             pdpc,eps0,om,fac1,fac2,fac3,pglo,ids,
     *             nzapjet,lzapjet)
          else if(mast(26).eq.3)then
               print*,' cyclt1:potential integrate along',
     *			  ' a geomagnetic field - 2nd var'

               call potnewmf2(delta,uts,pkp,bmpy,bmpz,
     *			 vsol,csol,ntr,ap,
     *             dst,al,ae,au,rads,nh,ddolgs,ddolgt,dtets,dtett,
     *             isp,readfl,kpars,kdf,kdu,nr,ldor,nfile,md,idt,
     *             its,nl2,par,pole,potef,mast,tet1,tet2,tet3,
     *             pdpc,eps0,om,fac1,fac2,fac3,pglo,ids,
     *             nzapjet,lzapjet)
          else
               print 234,mast(26)
234            format(' cyclt1 : ********  ERROR ! ********** '/
     *      '       mast(26)=',i3,' !!!!!  Stop !')

          stop
          endif
    8   continue
        nin=0
    3   if(utt.ge.uts) go to 4
          nf=6
          if(nin.ne.0) nf=13

            if(mass(13).ne.0) then
              call rvert(nf,kpart,ddolgt,ntsl,nl,kdf,ldor,isp,
     *                   pole,nv,idt,vert,mast)
              call cyclt2_bas(mast,ntsl,nl,par,nr,pari,ni,par1,PAR2,b,c,
     *                    its,park,ks,gins,rads,nh,ddolgt,dtett,ddolgs,
     *                    dtets,kpart,int,ins,rmaxt,ntr,nv,idt,ids,
     *                    dtt0,dtt,potef,q,u,nl2,vert,utt,kdu,kdf,
     *                    ldor,isp,pole,nin,verno,vdr,sole,solen,nse,
     *                    qom,qmax,iqo,mass)
            end if
            nin=1
            utt=utt+dtt
            if(utt.le.uts) go to 5
              utt=uts
    5       continue
            go to 3
    4   continue
        nin=2
        nf=13
        if(mass(13).ne.0) then
         call rvert(nf,kpart,ddolgt,ntsl,nl,kdf,ldor,isp,
     *               pole,nv,idt,vert,mast)
         call cyclt2_bas(mast,ntsl,nl,par,nr,pari,ni,par1,PAR2,b,c,its,
     *                park,ks,gins,rads,nh,ddolgt,dtett,ddolgs,dtets,
     *                kpart,int,ins,rmaxt,ntr,nv,idt,ids,dtt0,dtt,
     *                potef,q,u,nl2,vert,utt,kdu,kdf,ldor,isp,pole,
     *                nin,verno,vdr,sole,solen,nse,qom,qmax,iqo,mass)
        end if
        i=i+1
        go to 1
    2 continue
      dtt0=dtt
      if(dtt0.le.dts) go to 7
        dtt0=dts
    7 continue
      tx1=(ut1+0.5)/3600.
      ntx1=tx1
      tx2=(tx1-ntx1)*60.
      ntx2=tx2
      tx3=(tx2-ntx2)*60.
      ntx3=tx3
      print 666
 666  format(' cyclt1:   beginning of copying: 5 and 6 files  ')
      nfile=5
      readfl=.false.
      npgl=kpars*nh*its*ids
      nkd=20
      nf=4
      call readwr(readfl,nfile,kdf,kdu,nkd,npgl,pglo,ldor,isp,nf)

   !!!!!!!!!!!!!!!!!!!!
      if(mass(13).ne.0) then
         call copmd(13,6,kdf,isp,ldor,kdu,par,nr,mass,mast)
      end if
      t1=(ut1+tau+0.5)/3600.
      nt1=t1
      t2=(t1-nt1)*60.
      nt2=t2
      t3=(t2-nt2)*60.
      nt3=t3
      print 779,nt1,nt2,nt3
  779 format('  ',i2,'.',i2,'.',i2)
      readfl=.false.
      ut3=ut1+dts
      call zait(readfl,5,nzapt,nzaps,nadrt,nadrs,isp,
     *          ldor,kdf,god,day,ut0,ut3,pole)
      if(mass(13).ne.0) then
        call zait(readfl,6,nzapt,nzaps,nadrt,nadrs,isp,
     *           ldor,kdf,god,day,ut0,ut3,pole)
      end if
c . . . writing date and time in inf-record
      read(4,rec=1) pole
      pole(2)=god
      pole(3)=day
      pole(4)=ut0
      pole(5)=ut3
      write(4,rec=1) pole
c . . .
      print 667
 667  format(' cyclt1:   end of copying:       5 and 6 files ')
c-------------------------------------------------------------------
      if(mast(31).ne.0)then
        print *,' writing in ouf-files'
        do 300 n=1,koob
          namefl=naou(n)
          ip=iput(n)
          nrecl=keut(n)
          nut=kut(n)
          call intut(uts,tut,sut,dut,hut,n,nut,
     *               et,ed,eh,ieisc,isat)
          if(ieisc.ne.0) then
            izap(n)=izap(n)+1

            call eiscat(et,ed,eh,namefl,par,par1,potef,pole,kpars,
     *                kpart,nh,
     *                its,ntr,idt,nl2,ldor,ntsl,nl,kdf,kdu,isp,tau,
     *                ks,ddolgs,dtets,dtett,nzaps,
     *                nadrs,nr,uts,park,dts,izap(n),ip,nrecl,mast)
          end if
  300   continue
        open(9,file='fiza',status='old')
        rewind 9
        write(9,'(10i5)')izap
        close(9)
      end if
      
      print*,'GSMTIP: cyclt1 - end'
      print 670
670   format(' *****  END OF TIME STEP  *****')
      CALL CPU_TIME(TT2)
      print*,' CALCULATING step TIME = ', tt2-tt1
      deallocate (E0,FAE,vert,qom)
!      stop
      return
      end
