c    VER _04_2014 DINAMIC MASSIVES
C    version with Ohot
      subroutine cyclt2_bas(mast,ntsl,nl,par,nr,pari,ni,par1,PAR2,b,c,its,
     *                  park,ks,gins,rads,nh,ddolgt,dtett,ddolgs,
     *                  dtets,kpart,int,ins,rmaxt,ntr,nv,idt,ids,
     *                  dtt0,dtt,potef,q,u,nl2,vert,utt,kdu,kdf,
     *                  ldor,isp,pole,nin,verno,vdr,sole,solen,nse,
     *                  qom,qmax,iqo,mass)
      dimension ntsl(nl),par2(nr),par1(nr),par(nr),pari(ni),park(ks),
     *          qom(nl),rads(nh),potef(ntr,idt,nl2),q(nv,nl),u(nl),
     *          kdu(20),kdf(20),vdr(ks),sole(nse),solen(nse),
     *          pole(ldor/4),vert(kpart,nv),
     *          gins(ins,nh,its,ids),mass(30),mast(40)
      allocatable cio1(:),cih1(:),cihe1(:),vio1(:),vih1(:),
     *          vihe1(:),ti1(:),te1(:),cio(:),
     *          cih(:),cihe(:),vio(:),vih(:),vihe(:),
     *          ti(:),te(:),co2(:),cn2(:),co(:),
     *          ch(:),che(:),cim(:),tn(:),vnq(:),
     *          vnu(:),vnv(:),qo(:),ht(:),tt(:),vdv(:),
     *          cOhot(:),Tohot(:),
     *          vdu(:),qsm(:)

      integer verno
      logical readfl
      data pi/3.14159265/

      allocate( cio1(nv),cih1(nv),cihe1(nv),vio1(nv),vih1(nv),
     *          vihe1(nv),ti1(nv),te1(nv),cio(nv),
     *          cih(nv),cihe(nv),vio(nv),vih(nv),vihe(nv),
     *          ti(nv),te(nv),co2(nv),cn2(nv),co(nv),
     *          ch(nv),che(nv),cim(nv),tn(nv),vnq(nv),
     *          vnu(nv),vnv(nv),qo(nv),ht(nv),tt(nv),vdv(nv),
     *          cOhot(nv),Tohot(nv),
     *          vdu(nv),qsm(nv))
      readfl=.true.
      md=1
      dolm=0.
      print *,'cyclt2-begin'
      do  i=1,idt
        fi=dolm/180.*pi
        if(mast(13).eq.0) go to 4
        mi=1
        do j=1,nl
          ne=3
          m=2
          call select_oh(ne,j,m,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *                ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *                qo,qsm,ht,tt,vdv,vdu,park,ks,cOhot,Tohot,nv)
     		call vdrift(i,j,idt,ntsl,nl,ddolgt,dtett,fi,ht,tt,potef,
     *                vdu,vdv,ntr,nl2)
          ne=2
          m=2


          call fillin(ne,j,m,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *                ti1,te1,vdv,vdu,vdr,ks)
          mast24=mast(24)
c
          call dreifGSM(j,ntsl,nl,ht,tt,dolm,vdv,vdu,dtt,
     *               par,nr,par1,PAR2,
     *               park,ks,q,u,kpart,ddolgt,kdf,ldor,isp,pole,nv,
     *               rmaxt,vert,nh,ntr,rads,nin,mi,mast)
c
          mi=2
        end do


        go to 5
c
    4   continue
      nfile=13
      if(nin.eq.0) nfile=6
c********
      readfl=.true.
      nomsl=1
      md=1
!	PRINT*, KPART,DOLM,DDOLGT,NOMSL,LDOR,ISP
!	PAUSE
      call wwt(readfl,nfile,kpart,dolm,ddolgt,nomsl,ntsl,nl,kdf,ldor,
     *         isp,md,par1,pole,nr,mast)

!		 print*,' after wwt'
!		 print *,(par1(j),j=1,1120)
!		 print *,(par1(j),j=3361,4377)
!		 print*,' '
!		 print *,(par1(j),j=1121,2240)
!		 print *,(par1(j),j=4378,5314)
!		 print *,'  '
!		 print *,(par1(j),j=2241,3360)
!		 print *,(par1(j),j=5315,6187)
!
!		 pause
    5  continue
        nfile=8
        kpar=int
        call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
     *           kdf,ldor,isp,md,pari,pole,ni,mast)
        call trubka_oh(ntsl,nl,pari,ni,park,ks,vdr,par1,nr,mast,dolm,nv,
     *              ddolgt,kdf,ldor,isp,pole,kpart,cio1,cih1,cihe1,
     *              vio1,vih1,vihe1,ti1,te1,cio,cih,cihe,vio,vih,vihe,
     *              ti,te,co2,cn2,co,ch,che,tn,qo,qsm,cim,vnq,vnu,vnv,
     *              ht,tt,vdv,vdu,dtt,sole,solen,nse,int,qom,qmax,iqo
     *              ,utt,cOhot,Tohot,mass)

        if(nin.eq.2) then
          call ints(dolm,par1,nr,rads,nh,ni,pari,ins,its,park,ks,
     *              ntsl,nl,ntr,kdf,ldor,isp,pole,kpart,vdr,dtett,
     *              u,ddolgs,dtets,gins,ids)
        end if
        dolm=dolm+ddolgt
      end do
      call bongi(gins,ins,nh,its,ids)
      nfr=14
      nfw=13
      call copmd(nfr,nfw,kdf,isp,ldor,kdu,par,nr,mass,mast)
      print *,'cyclt2 - end'
      deallocate( cio1,cih1,cihe1,vio1,vih1,
     *          vihe1,ti1,te1,cio,
     *          cih,cihe,vio,vih,vihe,
     *          ti,te,co2,cn2,co,
     *          ch,che,cim,tn,vnq,
     *          vnu,vnv,qo,ht,tt,vdv,
     *          cOhot,Tohot,
     *          vdu,qsm)
      return
      end


