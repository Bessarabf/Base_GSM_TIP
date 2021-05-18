! subroutine cyclt1_bas, aepp, magsm, aeppVY, aurprecip, intut, rvert, 
!     inpout, eiscat, wws, potvol
!  function pefvol
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
      
      call shar_bas(day,god,dayt,godt,uts,tau,dts,sole,solu,solen,nsu,
     *           nse,kpars,ins,int,rads,nh,gkoor,its,ids,ddolgs,
     *           dtets,dtett,fa,fs,ap,pkp,dst,ae,al,au,bmpz,bmpy,
     *           vsol,ps,csol,mass,delta,kdu,kdf,ldor,isp,pole,par,
     *           nr,pari,par1,ni,ddolgt,ntsl,nl,verno,park,ks,potef,
     *           nl2,ntr,gins,solet,ut0,qom,qmax,iqo,mast,pglo,
     *           pril,kpa,nt,E0,FAE)
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
      deallocate (E0,FAE,vert,qom)
      return
      end
!------------------------------------------------------------------------------
	subroutine aepp(mass,Kp,nl,idt,ddolgs,dtet,ut,del0,E0,Q,its,ids)
     
      dimension Af(6),Bf(6),Cf(6),Df(6),Ae(6),Be(6),Ce(6),De(6)

      dimension E0(its,ids),Q(its,ids),mass(30)
     	
      dimension Acf(6,6),Bcf(6,6),Ccf(6,6),Dcf(6,6),Asf(6,6),Bsf(6,6),
     *Csf(6,6),Dsf(6,6),Ace(6,6),Bce(6,6),Cce(6,6),Dce(6,6),Ase(6,6),
     *Bse(6,6),Cse(6,6),Dse(6,6)
      dimension Akf(6),Bkf(6),Ckf(6),Dkf(6),Ake(6),Bke(6),Cke(6),Dke(6)
      real Kp,Kpm1,Kpm2,Kp_model(6)

	allocatable  qe0(:,:)
     	allocate (qe0(its,ids))

	data pi/3.14159265359/
      data Kp_model/.75,2.25,3.75,5.25,7.,9./
	
      open(1,file='Zhang_Paxton.dat')
      read(1,1)(Akf(i),i=1,6)
      do j=1,6
        read(1,1)(Acf(j,i),i=1,6)
        read(1,1)(Asf(j,i),i=1,6)
      end do
      read(1,1)(Bkf(i),i=1,6)
      do j=1,6
        read(1,1)(Bcf(j,i),i=1,6)
        read(1,1)(Bsf(j,i),i=1,6)
      end do
      read(1,1)(Ckf(i),i=1,6)
      do j=1,6
        read(1,1)(Ccf(j,i),i=1,6)
        read(1,1)(Csf(j,i),i=1,6)
      end do
      read(1,1)(Dkf(i),i=1,6)
      do j=1,6
        read(1,1)(Dcf(j,i),i=1,6)
        read(1,1)(Dsf(j,i),i=1,6)
      end do
      read(1,1)(Ake(i),i=1,6)
      do j=1,6
        read(1,1)(Ace(j,i),i=1,6)
        read(1,1)(Ase(j,i),i=1,6)
      end do
      read(1,1)(Bke(i),i=1,6)
      do j=1,6
        read(1,1)(Bce(j,i),i=1,6)
        read(1,1)(Bse(j,i),i=1,6)
      end do
      read(1,1)(Cke(i),i=1,6)
      do j=1,6
        read(1,1)(Cce(j,i),i=1,6)
        read(1,1)(Cse(j,i),i=1,6)
      end do
      read(1,1)(Dke(i),i=1,6)
      do j=1,6
        read(1,1)(Dce(j,i),i=1,6)
        read(1,1)(Dse(j,i),i=1,6)
      end do
      close(1)
    1 format(6f12.7)
      dpi=pi+pi
c     Определение интервала для Kp индекса и Kpm1 и Kpm2
      if(Kp.le.0.75)then
        Kpm1=0.75
        Kpm2=2.25
	end if
      if(Kp.gt.0.75.and.Kp.le.9.)then
	  do i=1,5
	    if(Kp_model(i).lt.Kp.and.Kp_model(i+1).ge.Kp)then
		  Kpm1=Kp_model(i)
		  Kpm2=Kp_model(i+1)
	    end if
	  end do
	end if
	if(Kp.gt.9.)then
	  Kpm1=7.00
	  Kpm2=9.00
	end if
c     Нахождение Hp индекса и Hpm1 и Hpm2
	if(Kp.le.5.)Hp=38.66*exp(0.1967*Kp)-33.99
	if(Kp.gt.5.)Hp=4.592*exp(0.4731*Kp)+20.47
	if(Kpm1.le.5.)Hpm1=38.66*exp(0.1967*Kpm1)-33.99
	if(Kpm1.gt.5.)Hpm1=4.592*exp(0.4731*Kpm1)+20.47
	if(Kpm2.le.5.)Hpm2=38.66*exp(0.1967*Kpm2)-33.99
	if(Kpm2.gt.5.)Hpm2=4.592*exp(0.4731*Kpm2)+20.47
	if(Kpm1.eq.0.75)then
        i=1
        ip=2 
      end if
	if(Kpm1.eq.2.25)then
        i=2
        ip=3 
      end if
	if(Kpm1.eq.3.75)then
        i=3
        ip=4 
      end if
	if(Kpm1.eq.5.25)then
        i=4
        ip=5 
      end if
	if(Kpm1.eq.7.00)then
        i=5
        ip=6 
      end if
c     Определение углов
        do j=1,nl
          tet=(j-1)*dtet
          alat=90.-tet
          x=90.-abs(alat)  
          do l=1,idt
            if(j.gt.10.and.j.lt.28)then
              E0(j,l)=1.e-3
		  qe0(j,l)=0.
              Q(j,l)=0.  
              goto 2 
            end if  
            dolg=(l-1)*ddolgs
            fi=dolg/180.*pi
            call magsm(ut,del0,fi,phism,1)
            tau=pi+phism
            if(tau.ge.dpi)tau=tau-dpi
            if(tau.lt.0.)tau=tau+dpi
c     Нахождение коэффициентов Af1-Df1 Af2-Df2 Ae1-De1 Ae2-De2
            do k=i,ip
              Af(k)=Akf(k)
              Bf(k)=Bkf(k)
              Cf(k)=Ckf(k)
              Df(k)=Dkf(k) 
              Ae(k)=Ake(k)
              Be(k)=Bke(k)
              Ce(k)=Cke(k)
              De(k)=Dke(k)
              do m=1,6
                Af(k)=Af(k)+Acf(k,m)*cos(m*tau)+Asf(k,m)*sin(m*tau) 
                Bf(k)=Bf(k)+Bcf(k,m)*cos(m*tau)+Bsf(k,m)*sin(m*tau) 
                Cf(k)=Cf(k)+Ccf(k,m)*cos(m*tau)+Csf(k,m)*sin(m*tau) 
                Df(k)=Df(k)+Dcf(k,m)*cos(m*tau)+Dsf(k,m)*sin(m*tau) 
                Ae(k)=Ae(k)+Ace(k,m)*cos(m*tau)+Ase(k,m)*sin(m*tau) 
                Be(k)=Be(k)+Bce(k,m)*cos(m*tau)+Bse(k,m)*sin(m*tau) 
                Ce(k)=Ce(k)+Cce(k,m)*cos(m*tau)+Cse(k,m)*sin(m*tau) 
                De(k)=De(k)+Dce(k,m)*cos(m*tau)+Dse(k,m)*sin(m*tau)
              end do 
            end do
c     Нахождение Eom1 E0m2 Qm1 Qm2
           E0m1=Ae(i)*exp((x-Be(i))/Ce(i))/(1.+exp((x-Be(i))/De(i)))**2
            E0m2=Ae(ip)*exp((x-Be(ip))/Ce(ip))/(1.+exp((x-Be(ip))/
     *        De(ip)))**2
           Qm1=Af(i)*exp((x-Bf(i))/Cf(i))/(1.+exp((x-Bf(i))/Df(i)))**2
            Qm2=Af(ip)*exp((x-Bf(ip))/Cf(ip))/(1.+exp((x-Bf(ip))/
     *        Df(ip)))**2
c     Нахождение E0 - энергии
		f1=(Kpm2-Kp)/(Kpm2-Kpm1)
		f2=(Kp-Kpm1)/(Kpm2-Kpm1)
		E0(j,l)=f1*E0m1+f2*E0m2

c     Нахождение Q - потока
		f1=(Hpm2-Hp)/(Hpm2-Hpm1)
		f2=(Hp-Hpm1)/(Hpm2-Hpm1)
		qe0(j,l)=f1*Qm1+f2*Qm2
            Q(j,l)=qe0(j,l)/E0(j,l)*6.24e8 
    2     continue
	  end do
	end do
        if(mass(30).eq.1)then	
          open(1,file='epres',access='direct',recl=723)
	    nrec=nrec+1
          write(1,rec=nrec)ut
	    nrec=nrec+1
          write(1,rec=nrec)Kp
	    nrec=nrec+1
          write(1,rec=nrec)Hp
	    do j=1,10
	      nrec=nrec+1
	      write(1,rec=nrec)(E0(j,l),l=1,24)
	    end do
	    do j=1,10
	      nrec=nrec+1
	      write(1,rec=nrec)(qe0(j,l),l=1,24)
	    end do
	    do j=1,10
	      nrec=nrec+1
	      write(1,rec=nrec)(Q(j,l),l=1,24)
	    end do
        end if
      close(1)
	deallocate (qe0)
	return 
	end
!------------------------------------------------------------------------------
      subroutine magsm(ut,del,phid,phism,j)
      data g10/30103.6/,g11/-2016.5/,h11/5682.6/,pi2/6.283185/
  900 format(' ',10g12.4)
      sq=g11**2+h11**2
      sqq=sqrt(sq)
      sq=sqrt(sq+g10**2)
      st0=sqq/sq
      ct0=g10/sq
      sl0=-h11/sqq
      cl0=-g11/sqq
      hour=ut/3600.
      al=0.2618*(hour-12.)
      sal=sin(al)
      cal=cos(al)
      ssd=sin(del)
      csd=cos(del)
      x1=csd*cal
      y1=-csd*sal
      z1=(x1*cl0+y1*sl0)*ct0-ssd*st0
      y1=y1*cl0-x1*sl0
      xmut=12.-3.8197*atan2(y1,z1)
      fi=-0.2618*xmut
      if(j.lt.0) go to 1
      phism=phid-fi-pi2/2.
      if(phism.lt.0) phism=phism+pi2
      if(phism.ge.pi2) phism=phism-pi2
      go to 2
    1 continue
      phid=phism+fi
      if(phid.ge.pi2) phid=phid-pi2
      if(phid.lt.0) phid=phid+pi2
    2 continue
c     print 900,ut,del,phid,phism,j
      return
      end
!------------------------------------------------------------------------------
      subroutine aeppVY(mass,nl,idt,ddolgs,dtet,AL,Dst,ut,del0,E0,Q)
      !      nl=its, idt=ids
      dimension E0(nl,idt),Q(nl,idt),MASS(30)
      allocatable qe0(:,:),ebpp1(:),ebpp2(:),ebpp3(:),ebpp4(:),
     *            ebpp5(:),gmltm(:)
      allocate(qe0(nl,idt+1),ebpp1(idt+1),ebpp2(idt+1),ebpp3(idt+1),
     *         ebpp4(idt+1),ebpp5(idt+1),gmltm(idt+1))      
      data pi/3.14159265359/
      dpi=pi+pi
      if(AL.ge.-5.)AL=-5.
      do l=1,idt+1
	 
        l00=l
	  if(l.eq.idt+1)l00=1 

        dolg=(l-1)*ddolgs
        fi=dolg/180.*pi
        call magsm(ut,del0,fi,phism,1)
        tau=pi+phism
        if(tau.ge.dpi)tau=tau-dpi
        if(tau.lt.0.)tau=tau+dpi
        gmlt=tau/pi*12.
        gmltm(l)=gmlt
          if(gmlt.eq.0.)then
            b1e=64.64+0.01*al+1.25*1.e-6*al*al+.02*dst
            b2e=66.69+0.009*al+7.78*1.e-7*al*al+.022*dst
            b5e=70.89-0.002269047619*al+0.00175*dst
            b6e=71.74-0.000994*al-0.0038*dst
            b2i=66.25158347+0.005469936384*al-5.430946029e-7*al*al
     *      +0.02255*dst+0.361
            b4s=67.48+0.007*al+1.36e-6*al*al+0.048*dst
          end if
          if(gmlt.gt.0..and.gmlt.le.3.)then
            b1e=64.77+0.012*al+1.2*1.e-6*al*al+.03*dst
            b2e=67.1+0.01*al+.58*1.e-6*al*al+.025*dst
            b5e=73.2+0.003*al+.0026*dst
            b6e=74.2+0.004*al+.01*dst
            b2i=68.7886781+0.0060622409195*al-2.7154730145e-7*al*al
     *      +0.03035*dst
c            b4s=(71.25629735+0.006298636364*al-2.916666667e-7*al*al
c     *      -0.3755-0.00545*dst+67.48+0.007*al+1.36e-6*al
c     *      *al+0.048*dst)*.5
            b4s=70.19100384905+0.0064990462755*AL+5.3379164115E-007*AL
     *          *AL+0.024404761905*Dst
          end if
          if(gmlt.gt.3..and.gmlt.le.6.)then
            b1e=64.90363095+0.01525*al+6.547619048*1.e-6*al*al
     *      +0.056*dst
c            b2e=68.29+0.014*al+7.187*1.e-6*al*al+0.033625*dst
            b2e=67.69140422+0.01406017316*AL+7.187229437e-006*AL*AL
     *          +0.03076785714*Dst+0.1494642857
            b5e=75.85+0.0061*al+0.0043*dst
            b6e=77.2+0.00746*al+0.012*dst
            b2i=70.06127273+0.006654545455*al+0.03815*dst+0.9035
            b4s=71.25629735+0.006298636364*al-2.916666667e-7*al*al
     *          +0.03125595238*Dst+0.7564880952
c     *      -0.3755-0.00545*dst
          end if
          if(gmlt.gt.6..and.gmlt.le.9.)then
            ecps=66.129+0.0142*al+6.871*1.e-6*al*al+.0524*dst
            pcps=72.6485+0.01603396486*al+8.97319837*1.e-6*al*al+
     *      0.023827*dst
            ebps=73.35+0.01145*al+3.292*1.e-6*al*al+0.02354761905*dst
            pbps=77.71435386+0.01015437229*al+2.79004329*1.e-6*al*al
     *      +0.01985714286*dst
            pllb=78.99303084381+0.01654489177*al+5.851731602*1.e-6*al
     *      *al+0.02108571429*dst
          end if
          if(gmlt.gt.9..and.gmlt.le.12.)then
            ecps=67.21+0.006470909091*al+0.0219*dst 
            pcps=74.03+0.00579197*al-3.40909*1.e-8*al*al+0.0293*dst 
            ebps=76.2+0.00806*al+0.0295*dst 
            pbps=78.11+0.007922727273*al+0.0186*dst
            pllb=80.53+0.01210984848*al+2.276515152*1.e-6*al*al+0.0293
     *      *dst
          end if
          if(gmlt.gt.12..and.gmlt.le.15.)then
            ecps=69.81670238+0.008624285714*al+2.071428571*1.e-6*al*al
     *      +0.04185*dst 
            pcps=74.2+0.01015458874*al+2.336580087*1.e-6*al*al
     *      +0.04185*dst 
            ebps=77.97+0.0137058807*al+3.324815877*1.e-6*al*al
     *      +0.0511*dst   
            pbps=79.984+0.01357380952*al+4.130952381*1.e-6*al*al
     *      +0.04585*dst
            pllb=81.18+0.01602409091*al+4.526515152*1.e-6*al*al
     *      +0.04745*dst
          end if
          if(gmlt.gt.15..and.gmlt.le.18.)then
            ecps=70.22+0.01153809524*al+5.345238095*1.e-6*al*al
     *      +0.06485*dst 
            pcps=73.266+0.01364214876*al+6.432113341*1.e-6*al*al
     *      +0.04532857143*dst 
            ebps=74.4+0.01529090909*al+6.735930736*1.e-6*al*al
     *      +0.0407*dst   
            pbps=78.+0.01693571429*al+8.94047619*1.e-6*al*al
     *      +0.02945*dst
            pllb=80+0.02479285714*al+1.272619048*1.e-5*al*al
     *      +0.02945*dst
          end if
          if(gmlt.gt.18..and.gmlt.le.21.)then
            b1e=68.91+0.012*al+3.85*1.e-6*al*al+0.05691428571*dst
c            b2e=71.14+0.011*al+3.61*1.e-6*al*al+0.05691428571*dst
            b2e=69.73307168+0.01133142191*AL+3.60955711E-006*AL*AL
     *          +0.9166517857+0.04044642857*Dst
            b5e=74.93+0.01275247934*al+8.820936639*1.e-6*al*al
     *      +0.0229*dst
            b6e=76.52+0.01780134525*al+1.155390033*1.e-5*al*al
     *      +0.01865*dst
            b2i=68.22612587+0.009239067599*al+2.610722611e-6*al*al
     *      +0.04798571429*dst+1.174142857
            b4s=71.48292992+0.01292863636*al+5.65530303e-6*al*al
     *      +0.024*dst+0.3
          end if
          if(gmlt.gt.21..and.gmlt.le.24.)then
            b1e=64.64+0.01*al+1.25*1.e-6*al*al+.02*dst
c            b2e=66.69+0.009*al+7.78*1.e-7*al*al+.022*dst
            b2e=66.71830828+0.009183808742*AL+7.78062795e-7*AL*AL
     *          +0.0226*Dst+0.07266666667
            b5e=70.89-0.002269047619*al+0.00175*dst
            b6e=71.74-0.000994*al-0.0038*dst
            b2i=66.25158347+0.005469936384*al-5.430946029e-7*al*al
     *      +0.02255*dst+0.361
c            b4s=67.48+0.007*al+1.36e-6*al*al+0.048*dst
            b4s=68.39+0.007*AL+1.36e-6*AL*AL+0.017*Dst
          end if
          if(gmlt.ge.0..and.gmlt.le.6..or.gmlt.gt.18..
     *    and.gmlt.le.24.)then
            ebpp1(l)=b1e
            ebpp2(l)=b2e
            ebpp3(l)=b4s
            ebpp4(l)=b5e
            ebpp5(l)=b6e
          end if
          if(gmlt.gt.6..and.gmlt.le.18.)then
            ebpp1(l)=ecps
            ebpp2(l)=pcps
            ebpp3(l)=ebps
            ebpp4(l)=pbps
            ebpp5(l)=pllb
          end if
        call aurprecip(gmlt,al,f1,f2,
     *         f4,f5,e1,e2,e4,e5,eaop,
     *         edaz,esdp,faop,fdaz,fsdp)
        do j=1,nl/2+1
          tet=(j-1)*dtet

          gmlat=90.-tet
          if(gmlt.ge.0..and.gmlt.le.6..or.gmlt.gt.18..
     *       and.gmlt.le.24.)then
            if(gmlat.gt.b1e.and.gmlat.le.b2e)then
              qe0(j,l)=f1
              E0(j,l00)=e1
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.b2e.and.gmlat.le.b4s)then
              qe0(j,l)=f2
              E0(j,l00)=e2
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.b4s.and.gmlat.le.b5e)then
              qe0(j,l)=f4
              E0(j,l00)=e4
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.b5e.and.gmlat.le.b6e)then
              qe0(j,l)=f5
              E0(j,l00)=e5
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
          end if
          if(gmlt.gt.6..and.gmlt.le.18.)then
            if(gmlat.gt.ecps.and.gmlat.le.pcps)then
              qe0(j,l)=fdaz
              E0(j,l00)=edaz
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.ge.ebps.and.gmlat.le.pbps)then
              qe0(j,l)=faop
              E0(j,l00)=eaop
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.pbps.and.gmlat.le.pllb)then
              qe0(j,l)=fsdp
              E0(j,l00)=esdp
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
          end if

        end do
      end do
! south hemisphere 
      
      do j=1,nl/2+1
        jc=nl-j+1
	  do l=1,idt+1
           l00=l
           if(l.eq.idt+1)l00=1 
           E0(jc,l00)=E0(j,l00)
           Q(jc,l00)=Q(j,l00)
        end do
      end do
		  
        if(mass(30).eq.1)then	
          open(1,file='epres',access='direct',recl=1371)
          nrec=nrec+1
          write(1,rec=nrec)ut
	      nrec=nrec+1
          write(1,rec=nrec)AL
	      nrec=nrec+1
          write(1,rec=nrec)Dst
	      do j=1,nl/2+1
	        nrec=nrec+1
	        write(1,rec=nrec)(E0(j,l),l=1,idt)
	      end do
	      do j=1,nl/2+1
	        nrec=nrec+1
	        write(1,rec=nrec)(qe0(j,l),l=1,idt)
	      end do
	      do j=1,nl/2+1
	        nrec=nrec+1
	        write(1,rec=nrec)(Q(j,l),l=1,idt)
	      end do
          close(1)
          open(1,file='bound',access='direct',recl=147)
	      nrec=nrec+1
          write(1,rec=nrec)ut
	      nrec=nrec+1
          write(1,rec=nrec)AL
	      nrec=nrec+1
          write(1,rec=nrec)Dst
	      do l=1,idt
	        nrec=nrec+1
            write(1,rec=nrec)gmltm(l),ebpp1(l),ebpp2(l),ebpp3(l),
     *          ebpp4(l),ebpp5(l)
	      end do
          close(1)
        end if

      print*,' '
      print 100,(ebpp1(i),i=1,idt)
      print*,' '
  100 format(' ',5f9.3)
      deallocate(qe0,ebpp1,ebpp2,ebpp3,ebpp4,ebpp5,gmltm)      
      return
      end
!------------------------------------------------------------------------------
      subroutine aurprecip(gmlt,al,f1,f2,
     *       f4,f5,e1,e2,e4,e5,eaop,edaz,esdp,
     *       faop,fdaz,fsdp)
      if(gmlt.eq.0.)then
        f1=exp(0.6571873362*alog(abs(al))-3.393664639) 
        f2=exp(0.709435388*alog(abs(al))-2.682839591)
        f4=exp(0.4970209835*alog(abs(al))-0.8709023689)
        f5=exp(0.1248129607*alog(abs(al))-1.84467514)
        e1=exp(0.2783876967*alog(abs(al))-0.553631512) 
        e2=exp(0.3491325617*alog(abs(al))-0.5661737632)
        e4=exp(0.2860529128*alog(abs(al))-0.6452854749)
        e5=exp(0.377738*alog(abs(al))-2.183403144)
      end if
      if(gmlt.gt.0..and.gmlt.le.3.)then
        f1=exp(0.5010929633*alog(abs(al))-2.057026824)
        f2=exp(0.7668563524*alog(abs(al))-2.961112584)
        f4=exp(0.5390297179*alog(abs(al))-1.599382872)
        f5=4.254511278E-5*abs(AL)+0.2174526316
c        f5=abs(f5) !!! correction 28.08.2015
        faop=exp(0.7229942079*alog(abs(AL))-2.72199829)
        e1=exp(0.208705056*alog(abs(al))+0.1701497954)
        e2=exp(0.2827463961*alog(abs(al))-0.2616186626)
        e4=exp(0.2847802859*alog(abs(al))-1.017263627)
        e5=exp(0.2500022797*alog(abs(al))-1.712236379)
        eaop=0.8683157161*alog(abs(AL))-2.02028153
      end if
      if(gmlt.gt.3..and.gmlt.le.6.)then
c        f1=exp(0.5010929633*alog(abs(al))-2.057026824)
c        f2=exp(0.7668563524*alog(abs(al))-2.961112584)
c        f4=exp(0.5390297179*alog(abs(al))-1.599382872)
c        f5=4.254511278e-5*al+0.2174526316
c        f5=abs(f5) !!! correction 28.08.2015
        f2=exp(0.8159742601*alog(abs(AL))-3.195652962)
        f4=exp(0.7038019846*alog(abs(AL))-3.54303041)
        f1=exp(0.4158926226*alog(abs(AL))-1.268071995)
        f5=exp(0.1468195661*alog(abs(AL))-2.846272624)
        faop=exp(0.8862633552*alog(abs(AL))-3.878830031)
        e1=exp(0.1667090503*alog(abs(al))+0.6541825456)
        e2=exp(0.198093427*alog(abs(al))+0.1380986798)
        e4=exp(0.2812593182*alog(abs(al))-1.615408825)
        e5=exp(-0.0523166127*alog(abs(al))-0.3659097527)
        eaop=0.7665735864*alog(abs(AL))-1.723066197
      end if
      if(gmlt.gt.6..and.gmlt.le.9.)then
        eaop=exp(0.08592111476*alog(abs(al))-0.2500901102)
        edaz=exp(0.1915357972*alog(abs(al))+0.5662103641)
        esdp=-5.e-5*abs(al)+0.3158333333
        faop=exp(0.4138865686*alog(abs(al))-1.567065282)
        fdaz=exp(0.3313262308*alog(abs(al))-1.443740347)
        fsdp=0.0007166666667*abs(al)+0.4697222222
c        fsdp=abs(fsdp) !!! correction 28.08.2015
      end if
      if(gmlt.gt.9..and.gmlt.le.12.)then
        eaop=exp(0.08496450811*alog(abs(al))-0.5313937638) 
        edaz=exp(0.2007310553*alog(abs(al))+0.7161005602)
        esdp=-1.e-5*abs(AL)+0.2756111111
        faop=exp(0.1615595083*alog(abs(al))-0.7343816142)
        fdaz=exp(0.07503442233*alog(abs(al))-0.7564862408)
c        fsdp=0.0006016666667*al+0.5025833333
c        fsdp=abs(fsdp) !!! correction 28.08.2015
cc       fsdp0=0.2638098631*alog(abs(al))-1.826307057
        fsdp=0.1815317224*alog(abs(AL))-0.2868825239
      end if
      if(gmlt.gt.12..and.gmlt.le.15.)then
        eaop=exp(0.06244459649*alog(abs(al))-0.8492598268) 
        edaz=exp(0.04497766381*alog(abs(al))+1.217024722)  
        esdp=5.142857143e-5*abs(al)+0.2312380952
        faop=exp(0.3052293503*alog(abs(al))-1.577706858)  
        fdaz=exp(-0.145536954*alog(abs(al))-1.148174158)
c        fsdp=0.00082*al+0.409
c        fsdp=abs(fsdp) !!! correction 28.08.2015
        fsdp=exp(0.3548040564*alog(abs(AL))-2.341480608)
      end if
      if(gmlt.gt.15..and.gmlt.le.18.)then
        eaop=exp(0.2654109679*alog(abs(al))-1.351216348) 
        edaz=exp(-0.05754892011*alog(abs(al))+1.149209589)  
        esdp=-2.285714286e-5*abs(al)+0.3101904762
        faop=exp(0.5976870883*alog(abs(al))-2.775956307)  
        fdaz=exp(-0.1464972011*alog(abs(al))-1.224727104) 
        fsdp=0.001574285714*abs(al)+0.3127142857
c        fsdp=abs(fsdp)  !!! correction 28.08.2015
      end if
      if(gmlt.gt.18..and.gmlt.le.21.)then
        f1=exp(0.5024018313*alog(abs(al))-3.406987641)
        f2=exp(0.5576887485*alog(abs(al))-3.099012084)
        f4=exp(0.7512055891*alog(abs(al))-2.799803464)
        f5=exp(0.4294515713*alog(abs(al))- 3.602047784)
        faop=exp(0.806511179*alog(abs(AL))-3.313109732)
        e1=exp(0.1264477554*alog(abs(al))-0.265670415)
        e2=exp(0.1617507433*alog(abs(al))-0.28674081)
        e4=exp(0.4355726882*alog(abs(al))-1.732201529)  
        e5=exp(0.5313973653*alog(abs(al))-3.340605919)
        eaop=0.385046894*alog(abs(AL))-1.365733462
      end if
      if(gmlt.gt.21..and.gmlt.le.24.)then
        f1=exp(0.6571873362*alog(abs(al))-3.393664639) 
        f2=exp(0.709435388*alog(abs(al))-2.682839591)
        f4=exp(0.4970209835*alog(abs(al))-0.8709023689)
        f5=exp(0.1248129607*alog(abs(al))-1.84467514)
        faop=3.504970748*alog(abs(AL))-13.14906522
        e1=exp(0.2783876967*alog(abs(al))-0.553631512) 
        e2=exp(0.3491325617*alog(abs(al))-0.5661737632)
        e4=exp(0.2860529128*alog(abs(al))-0.6452854749)
        e5=exp(0.377738*alog(abs(al))-2.183403144)
        eaop=1.029401751*alog(abs(AL))-2.632421127
      end if
      return
      end   
!------------------------------------------------------------------------------
       subroutine intut(uts,tut,sut,dut,hut,n,nut,
     *                  et,ed,eh,ieisc,isat)

      dimension tut(isat,100),sut(isat,100),dut(isat,100),hut(isat,100)
      ut=uts/3600.
      if(ut.lt.tut(1,n).and.(ut+24.).gt.tut(nut,n))then
        ieisc=0
        return
      end if
      i=1
    1 continue
      ip=i+1
      if(ip.gt.nut)then
        ieisc=0
      else
        if(ut.ge.tut(i,n).and.ut.le.tut(ip,n))then
          ieisc=1
          if(nut.eq.2)then
            et=sut(1,n)
            ed=dut(1,n)
            eh=hut(1,n)
          else
            a=tut(i,n)
            b=tut(ip,n)
            del=(ut-a)/(b-a)
            a=sut(i,n)
            b=sut(ip,n)
            et=a+(b-a)*del
            a=dut(i,n)
            b=dut(ip,n)
            ed=a+(b-a)*del
            a=hut(i,n)
            b=hut(ip,n)
            eh=a+(b-a)*del
          end if
        else
          if((ut+24.).ge.tut(i,n).and.(ut+24.).le.tut(ip,n))then
            ieisc=1
            if(nut.eq.2)then
              et=sut(1,n)
              ed=dut(1,n)
              eh=hut(1,n)
            else
              a=tut(i,n)
              b=tut(ip,n)
              del=(ut+24.-a)/(b-a)
              a=sut(i,n)
              b=sut(ip,n)
              et=a+(b-a)*del
              a=dut(i,n)
              b=dut(ip,n)
              ed=a+(b-a)*del
              a=hut(i,n)
              b=hut(ip,n)
              eh=a+(b-a)*del
            end if
          else
            i=ip
            goto1
          end if
        end if
      end if
      return
      end
!------------------------------------------------------------------------------
      subroutine rvert(nf,kpart,ddolgt,ntsl,nl,kdf,ldor,isp,
     *       pole,nv,idt,vert,mast)
      dimension pole(ldor/4),kdf(20),vert(kpart,nv)
     *       ,ntsl(nl),mast(40)
      allocatable pp(:,:)

	logical readfl
      
	allocate (pp(kpart,nv))

	npp=kpart*nv
      readfl=.true.
      md=0
      nomsl=1
      dolm=0.
 !     do 8 i=1,nv
 !       do 9 j=1,kpart
          vert=0.
 !   9   continue
 !   8 continue
      do 10 k=1,idt
        call wwt(readfl,nf,kpart,dolm,ddolgt,nomsl,ntsl,nl,
     *           kdf,ldor,isp,md,pp,pole,npp,mast)
        do 11 i=1,nv
          do 12 j=1,kpart
            vert(j,i)=vert(j,i)+pp(j,i)
   12     continue
   11   continue
        dolm=dolm+ddolgt
   10 continue
      do 13 i=1,nv
        do 14 j=1,kpart
          vert(j,i)=vert(j,i)/float(idt)
   14   continue
   13 continue
      deallocate (pp)
      return
      end
!------------------------------------------------------------------------------
!
! чтение - запись параметров трубки
! readfl - условие чтения - записи
! kdf - массив, определяющий начальную запись для чтения
! ldor - длина записи (байт)
! isp - номер записи
! nob - 1
! nf  - номер файла в open(nf
! nfile - номер элемента в kdf - начало записи shar, pot и т.д.
! lpar = nh*kpar 
! par - одномерный массив с параметрами, который заполняется
! pole - рабочий массив, сюда читаем данные
! nr -  
      subroutine inpout(readfl,nfile,kdf,ldor,
     *       isp,nob,lpar,par,pole,nr)
      logical readfl
      dimension par(nr),pole(ldor/4)
      dimension kdf(20)
!	if (nfile.eq.5) then
!	 print*,kdf
!	 stop
!	end if
      mdor=ldor/4
      nf=4
      if(nfile.gt.8) then 
         nf=5 ! файл f5
      !print*, ' inpout file5', readfl
	end if
!    integer number of records
      ndor=nob/mdor
!    difference of addresses
      nost=nob-ndor*mdor
!    
      nsvob=mdor-nost
      if(nost.eq.0)nsvob=0
      if(lpar.le.nsvob) then
        kd=0
        ln=0
        kz=0
        lost=0
      else
        ln=lpar-nsvob
        kd=ln/mdor
        kz=kd*mdor
        lost=ln-kz
      end if
      isp=ndor+kdf(nfile)+1
      if(readfl) go to 7 ! идти читать, иначе запись
! 
        if(nsvob.ne.0) then 
          read(nf,rec=isp) pole
          n=nost+1
          k=n+lpar-1
          if(nsvob.le.lpar)k=mdor
          i=1
          do j=n,k
            pole(j)=par(i)
            i=i+1
          end do 
          write(nf,rec=isp) pole 
        end if 
        if(kd.ne.0) then 
          n=nsvob+1
          k=n+mdor-1
          is=isp
          if(nsvob.ne.0)is=is+1
          do m=1,kd
            write(nf,rec=is)(par(i),i=n,k)
            n=n+mdor
            k=k+mdor
            is=is+1
          end do 
        end if 
        if(lost.ne.0) then
          is=isp
          if(nsvob.ne.0)is=isp+1
          if(kd.ne.0)is=is+kd
!          print*,'inpouT',' rec=',isp
          read(nf,rec=is) pole
          i=1
          n=nsvob+kz+1
          k=n+lost-1
          do j=n,k
            pole(i)=par(j)
            i=i+1
          end do
          write(nf,rec=is) pole
        end if 
        go to 6
    7 continue ! чтение
        if(nsvob.ne.0) then
          read(nf,rec=isp) pole
          isp=isp+1
          k=lpar
          if( nsvob.le.lpar) k=nsvob
          i=nost+1
          do j=1,k
            par(j)=pole(i)
            i=i+1
          end do
        end if 
        if(ln.ne.0) then
          n=nsvob+1
          k=ln+nsvob
          k1=n+mdor-1
          if(ln.le.mdor) k1=k
          m=kd
          if(lost.ne.0)m=m+1
          do j=1,m
            read( nf,rec=isp)(par(i),i=n,k1)
            isp=isp+1
            n=n+mdor
            k1=k1+mdor
            if(k1.gt.k)k1=k
          end do 
        end if
    6 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine eiscat(etet,edol,eheig,namefl,par1,par2,potef,pole,
     *                  kpars,kpart,nh,its,ntr,
     *                  idt,nl2,mdor,ntsl,nl,kdf,kdu,isp,tau,ks,
     *                  ddolg,dtets,dtet,nzap,adres,nr,ut,park,dts,
     *                  nrec,ip,nrecl,mast)

	parameter(kpars16=16)
	!kpars16 - only 16 parameters
      character*25 namefl
      logical readfl/.true./
      integer ntsl(nl),kdf(20),adres,kdu(20)
      
	real hm(70),cm1(70,2,8)
     *    ,fil(8,70)
      
      integer npar(8)/1,2,3,4,5,6,7,8/
      dimension par1(kpars,nh,its), par2(kpars,nh,its)
     *         ,pole(mdor/4), potef(ntr,idt,nl2),park(ks)
     *         ,nf(17),mast(40)
     
      allocatable a(:),cmi(:,:),ssh(:), out(:,:)

	
	common /flstat/ nf
      data  imin/1/, ns/1/
      data md/1/
	
      allocate (a(its),cmi(its,nh),ssh(its),out(its,nh))

      ntpl=ks/2
      dolg1=ifix(edol/ddolg)*ddolg
      if(dolg1.eq.360.)dolg1=0.
      dolg2=dolg1+ddolg
      if(dolg2.eq.360.)dolg2=0.
      open (12,file=namefl,access='direct',recl=nrecl)
          nfile = 5
      call wws (readfl,nfile,kpars,dolg1,ddolg,
     *          tet,dtets,nh,kdf,mdor,isp,md,par1,pole,nr)
      call wws (readfl,nfile,kpars,dolg2,ddolg,
     *          tet,dtets,nh,kdf,mdor,isp,md,par2,pole,nr)
          coef2=(edol-dolg1)/ddolg
          coef1=1.-coef2
      do 1 k=1,its
        do 1 j=1,nh
         do 10 i=1,6
          if(par1(i,j,k).lt.1.e-3) par1(i,j,k)=1.e-3
          if(par2(i,j,k).lt.1.e-3) par2(i,j,k)=1.e-3
          par1(i,j,k)=alog10(par1(i,j,k))
          par2(i,j,k)=alog10(par2(i,j,k))
 10      continue
         do 11 i=13,16
          if(par1(i,j,k).lt.1.e-3) par1(i,j,k)=1.e-3
          if(par2(i,j,k).lt.1.e-3) par2(i,j,k)=1.e-3
          par1(i,j,k)=alog10(par1(i,j,k))
          par2(i,j,k)=alog10(par2(i,j,k))
 11      continue
         do 1 i=1,kpars16
            par1(i,j,k)=par1(i,j,k)*coef1+par2(i,j,k)*coef2
 1      continue
          do 2 i=1,its
            a(i)=90.-(i-1)*dtets
        if ( etet.ge.a(i) ) then
          i1=i
          i2=i-1
          coef2=( etet-a(i1) )/dtets
          coef1=1.-coef2
          goto 3
        endif
 2    continue
 3    continue
        do 4 n=1,nh
          do 5 k=1,kpars16
            out(k,n)=par1(k,n,i1)*coef1+par1(k,n,i2)*coef2
            if(k.le.6.or.k.ge.13) out(k,n)=10.**(out(k,n))
 5        continue
 4      continue
          d0=edol
          tet=90.-etet
      idm=d0/ddolg
      coef1=(d0-idm*ddolg)/ddolg
      coef2=1.-coef1
      do 7 i=1,2
        id=idm+i-1
        if(id.ge.idt)id=id-idt
        dolg=ddolg*id
          nfile = 6
        call wwt(.true.,nfile,kpart,dolg,ddolg,i,ntsl,nl,kdf,mdor,
     *           isp,1,
     *           par1,pole,nr,mast)
        do 7 i1=1,kpart
!          call plosk(par1,park ,hm,cmi,ssh,tet,cm1(1,i,i1),ns,
!     *               ip,imin,ntsl,ntpl,nl,npar(i1),its,nh,nl2)
 7      continue
      do 8 i1=1,kpart
        do 8 id=1,ip
          fil(i1,id)=cm1(id,1,i1)*coef2+cm1(id,2,i1)*coef1
 8    continue
          do 9 j=1,ip
          do 9 i=1,3
            fil(i,j)=10**fil(i,j)
 9        continue
          nfile = nf(15)
        n3=ntr*idt*nl2
        call wpotef(.true.,potef,n3,kdf,kdu,mdor,isp)
        id1=idm+1
        if(id1.eq.idt+1)id1=1
        idd1=id1+1
        if(idd1.gt.idt)idd1=idd1-idt
        it1=tet/dtet+1
          rip = ip
        write(12,rec=nrec) ut,etet,edol,eheig,rip,
     *  ((OUT(I,J),I=1,KPARS16),J=1,NH),((FIL(I,J),I=1,KPART),J=1,Ip),
     *  potef(ntr,id1,it1),potef(ntr,id1,it1+1),potef(ntr,idd1,it1)
     *  ,potef(ntr,idd1,it1+1)
        close (12)
!        write(*,'(a8,a21,a8)')' end of ',namefl,' filling'
	deallocate (a,cmi,ssh,out)
      return
      end
!------------------------------------------------------------------------------
      subroutine wws (readfl,nfile,kpar,dolg,ddolgs,
     *       tet,dtets,nh,kdf,ldor,isp,md,par,pole,nr)
      dimension kdf(20),par(nr), pole(ldor/4)
      logical readfl
c     print 100
  100 format(' wws')
c     call timen
      nn=180./dtets+1
      lpar=nh*kpar
      nob1=0
      nob2=0
      if(dolg.ne.0.)nob2=lpar*nn*(dolg/ddolgs)
!
      if(md.eq.0) then 
        if(tet .ne.0.)nob1=(tet/dtets)*lpar
      else 
        lpar=lpar*nn
      end if
      nob=nob1+nob2
      call inpout(readfl,nfile,kdf,ldor,
     *       isp,nob,lpar,par,pole,nr)
c     call timen
      return
      end
!------------------------------------------------------------------------------
c . . . ver. 2012 
      subroutine potvol(delta,ut,kp,bmpy,alt,ddolg,dtet,potef,
     *                 idt,nl2,ntr,b,c,isp,readfl,kdf,kdu,ldor)
      double precision pi
      real kp
      dimension potef(ntr,idt,nl2),kdf(20),kdu(20)
      logical readfl
      data pi/3.14159265359d0/
      i2=idt/2
      print 7
      j1=(nl2+1)/2
      do 6 i=1,idt
        fi=(i-1)*pi/i2
        do 6 j=1,nl2
          if(j.ne.1)go to 1
            tet=0.
            go to 5
    1     continue
          if(j.ne.nl2)go to 2
            tet=pi
            go to 5
    2     continue
          if(j.ne.j1)go to 3
            tet=pi*0.5
            go to 5
    3     continue
          if(j.gt.j1)go to 4
            tet=(c-(j-j1+1)*dtet)*pi/180.
            go to 5
    4     continue
            tet=(b-(j-2)*dtet)*pi/180.
    5  continue
        potef(16,i,j)=pefvol(delta,ut,kp,bmpy,alt,tet,fi,j,j1,nl2)
    6 continue
    7 format(' potvol')
      na=16
      readfl=.false.
      nanl=na*idt*nl2
      call wpotef(readfl,potef,nanl,kdf,kdu,ldor,isp)
      
      pmaxn=potef(na,1,1)
      pminn=potef(na,1,1)
      pmaxs=potef(na,1,nl2)
      pmins=potef(na,1,nl2)
      do8j=1,j1
        k=nl2-j+1
        do8i=1,idt
          pn=potef(na,i,j)
          ps=potef(na,i,k)
          if(pmaxn.lt.pn)pmaxn=pn
          if(pminn.gt.pn)pminn=pn
          if(pmaxs.lt.ps)pmaxs=ps
          if(pmins.gt.ps)pmins=ps
    8 continue
      dpn=(pmaxn-pminn)*1.e-11
      dps=(pmaxs-pmins)*1.e-11
      print  9,dps,dpn
    9 format(' ',35x,'dps=',1pe9.2,' kV',10x,'dpn=',1pe9.2,' kV')
      return
      end
!------------------------------------------------------------------------------
       function pefvol(delta,utt,pkp,bmpy,alt,tet,fi,j,j1,nl2)
      real l,l0,l1,l3
      double precision pi
c     data fic1/8.6e11/,fic2/3.2e12/,fip/9.e12/,pi/3.14159265359d0/,
      data fic1/8.6e11/,fic2/4.0e12/,fip/9.e12/,pi/3.14159265359d0/,
     *l0/14.9/,l1/8.5/,l3/33.2/,re/6371.02e5/,
     *c1/1.908395e-1/,c2/4.6479083e-1/,c3/9.184622e-1/,
     *c4/2.029585/,c5/4.4961194e-2/
      pp=pi+pi
c     print 113
  113 format(' pefvol')
      if(j.eq.1.or.j.eq.nl2)go to 3
        l=(re+alt)/(re*sin(tet)**2)
c       tau = fi-fih(delta,utt)+pi
        call magsm(utt,delta,fi,phism,1)
        tau=pi+phism
        if(tau.ge.pp)tau=tau-pp
        if(tau.lt.0.)tau=pp+tau
        if(l.gt.l0)go to 3
          if(pkp.gt.2..and.pkp.le.4.)goto1
            pefvol=fic1*(l/l0)**2*sin(tau)
            return
    1     continue
          if(l.gt.l1)goto2
            tau=tau-c1
            pefvol=fic2*c2*(l/l1)**2*sin(tau)
            return
    2     continue
            tau=tau-c3+.34*alog(l)
            pefvol=fic2*(l/l0)**1.365*sin(tau)
            return
    3 continue
        if(bmpy.ne.0.)goto4
          pefvol=0.
          if(j.eq.1.or.j.eq.nl2)return
          c=0.
          go to 5
    4   continue
          c=fip
          if(bmpy.gt.0..and.j.gt.j1.or.bmpy.lt.0..and.j.lt.j1)c=-fip
          if(j.ne.1.and.j.ne.nl2)go to 5
            pefvol=c*c5
            return
    5     continue
          d=sqrt(l0/l)*sin(tau)
          if(pkp.gt.2..and.pkp.le.4.)goto6
            d=fic1*d
            goto7
    6     continue
            d=fic2*d
    7     continue
          if(l.gt.l3)goto8
            pefvol=c*c4*(sqrt(l/l0)-1.)**2/l+d
            return
    8     continue
            pefvol=c*(c5-1./l)+d
            return
      end