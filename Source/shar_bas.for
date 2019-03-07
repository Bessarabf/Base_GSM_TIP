c   shar_HAM - integrate with HAMMONIA 
c   ver.    18.05.18 Ion Drag sent to HAMMONIA 
c   ver.    10.05.18 Joul heating sent to HAMMONIA
      subroutine shar_bas(day,god,dayt,godt,uts,tau,dts,sole,solu,
     *           solen,nsu,nse,kpars,ins,int,rads,nh,gkoor,
     *           its,ids,ddolgs,dtets,dtett,fa,fs,ap,pkp,dst,
     *           ae,al,au,bmpz,bmpy,vsol,ps,csol,mass,delta,
     *           kdu,kdf,ldor,isp,pole,par,nr,pari,par1,ni,
     *           ddolgt,ntsl,nl,verno,park,ks,potef,nl2,ntr,
     *           gins,solet,ut0,qom,qmax,iqo,mast,pgl,pril,
     *           kpa,nt,E0,FAE)
!
!      USE mo_ham_gsm
      logical readfl
      integer god,day,godt,dayt,verno
    
      dimension sole(nse),solu(nsu),solet(nse),rads(nh),park(ks)
     *         ,gkoor(2,its,ids),kdu(20),kdf(20),pole(ldor/4)
     *         ,par(kpars,nh,its),pari(ins,nh,its),ps(10),qom(nl)
     *         ,potef(ntr,ids,nl2),par1(kpars,nh,its)
     *         ,pgl(kpars,nh,its,ids),solen(nse),ntsl(nl)
     *         ,gins(ins,nh,its,ids),mast(40),mass(30),
     *          pril(kpa,its,ids,nt)
      dimension E0(its,ids),FAE(its,ids)
      allocatable vir(:,:,:),vid(:,:,:)
     *           ,vim(:,:,:),parj(:,:,:),parj1(:,:)

c
      data nj/0/
      nj=0
      print 106
106   format (' GSMTIP: shar - start             ')
      
	allocate (vir(nh,its,ids),vid(nh,its,ids),
     *         vim(nh,its,ids),parj(nh,its,ids),parj1(nh,its))


      if(nj.eq.0) then
        call globr  (pgl,par,kpars,ddolgs,dtets
     *              ,nh,kdf,ldor,isp,pole,nr,its,ids,mass)
      end if
      if(mass(9).gt.3) mass(9)=mass(9)-3
      call zamion (gins,ins,its,ids,nh,uts,mass)
      do 10 k = 1 , nh
       do 10 i = 1 , its
        do 10 j = 1 , ids
         vir(k,i,j)=0.
         vim(k,i,j)=0.
         vid(k,i,j)=0.
  10  continue

      if(mass(3).eq.0) go to 11
c        if(mast(26).eq.2.or.mast(26).eq.3)then
c          call vmionmf(pgl,kpars,nh,its,rads,dtett,dtets,ddolgs
c     *          ,potef,ntr,nl2,ids,vim,vid,vir,mast,mass)
c        else
          call vmion (pgl,kpars,nh,its,rads,dtett,dtets,ddolgs
     *               ,potef,ntr,nl2,ids,vim,vid,vir,mast,mass)
c        end if
  11  continue
      
      call terpot_bas(day,god,dayt,godt,uts,tau,dts,solet,sole,solu
     *           ,nsu,nse,kpars,rads,nh,gkoor,its,ddolgs
     *           ,dtets,fa,fs,ap,pkp,dst,ae,al,au,bmpz,bmpy,
     *            mass,delta,pgl,gins,ids,ins,isp,
     *            vir,vid,vim,verno,parj,potef,ntr,nl2,pril,kpa,nt)
      PRINT*,' TERMOS END'
      
       mass(9)=mass(9)+1

c
c     call wrpot(potef,ntr,ids,nl2)
c
c     . . . atomic H calculated
      if(mass(13).ne.0.and.mast(32).ne.0) then
         print*,'atomh'
         call atomh(pgl,kpars,nh,its,ids,gkoor,rads,
     *              fs,fa,ap,uts,day)
      end if
      PRINT*, 'GSMTIP: glqint start'
      call glqint(pgl,par,pari,kpars,nh,its,ids,kdf,ldor,pole,
     *            isp,ddolgs,dtets,sole,solen,rads,gkoor,
     *            uts,delta,nse,gins,ins,mass,dts,ps,park,
     *            nr,ni,kpart,ks,ddolgt,int,ntsl,nl,parj,parj1,
     *            ntr,qom,qmax,iqo,mast,vir,E0,FAE)

      nj=1
      print 107
107   format(' GSMTIP: shar - finish           ')
      deallocate (vir,vid,vim,parj,parj1)
       return
       end
