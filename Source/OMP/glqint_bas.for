! Base version glqint 11/03-2019
!
      SUBROUTINE GLQINT_BAS(pgl1,par,pari,kpars,nh,its,ids,kdf,ldor,
     &                  pole,isp,ddolgs,dtets,sole,solen,rads,gkoor,
     &                  uts,delta,nse,gins,ins,mass,dts,ps,park,
     &                  nr,ni,kpart,ks,ddolgt,int_,ntsl,nl,parj,parj1,
     &                  ntr,qom,qmax,iqo,mast,vir,E0,FAE)
#ifdef _OPENMP
      USE OMP_LIB
#endif
      dimension pgl1(kpars,nh,its,ids),gins(ins,nh,its,ids)
     &         ,par(kpars,nh,its),pari(ins,nh,its)
     &         ,sole(nse),solen(nse),rads(nh),gkoor(2,its,ids)
     &         ,parj(nh,its,ids),parj1(nh,its),ps(10)
     &         ,mass(30),kdf(20),pole(ldor/4),ntsl(nl),park(ks)
     &         ,qom(nl),mast(40)
      DIMENSION E0(its,ids),FAE(its,ids)
      DIMENSION vir(nh,its,ids)
      LOGICAL readfl
      
      DATA key/0/

      readfl = .FALSE.
      IF ( Mass(20).EQ.1 ) CALL MOLIO2(Pgl1,Gins,Rads,Vir,Kpars,Nh,Its,
     &                                 Dts,Ids,Ddolgs,Ins,Dtets,Ntr)
      PRINT * , ' glqint '
!------------------------------------------------------------------------------
!!! longitude cycle
#ifdef _OPENMP
      CALL omp_set_num_threads(8)
#endif
!$OMP PARALLEL DEFAULT( SHARED )
!$OMP& FIRSTPRIVATE ( Ddolgs, Ddolgt, Delta, dolg, Dtets,
!$OMP& Dts, E0, Fae, Gkoor,  Ids, Ins, Int_, Iqo, Isp, Its, j, Kdf, key, Kpars,
!$OMP& Ks, Ldor, Mass, Mast, nfile, Nh, Ni, Nse, Nl, Nr, Ntr, Ntsl, Par, Pari, 
!$OMP& Parj1, Park, Pole, Ps, Qmax, Qom, Rads, Sole, Solen, Uts )
!$OMP DO SCHEDULE (DYNAMIC, 1)
      DO j = 1 , Ids
        dolg = Ddolgs*(j-1)
!!!form 2d massiv (for compatibility with old version GSM TIP)
        Par(1:Kpars,1:Nh,1:Its) = Pgl1(1:Kpars,1:Nh,1:Its,j)
        Pari(1:Ins,1:Nh,1:Its) = Gins(1:Ins,1:Nh,1:Its,j)
        CALL IONIZ_AE(Sole,Solen,Rads,Par,Parj1,Gkoor,Uts,Dtets,dolg,
     &                Ddolgs,Delta,Nh,Its,Ids,Nse,Kpars,Mass,Ps,E0,Fae)
!------------------------------------------------------------------------------
        IF ( Mass(20).EQ.0 ) THEN
          CALL MOLION(Par,Kpars,Nh,Its,Pari,Ins,Mass,Dts,Ntr)
        ELSEIF ( Mass(20).EQ.2 ) THEN
!	   call molio3S(cO2plus,cNOplus,par,ids,its,nh,kpars,pari,ins,
!     *                mass,dts,j,ntr,key)
        CALL MOLIO3NIT_BAS(Par,Ids,Its,Nh,Kpars,Pari,Ins,Mass,Dts,j,
     &                         Ntr,key)
        ENDIF ! Mass(20).EQ.0
!------------------------------------------------------------------------------
        PRINT *,'j=',j,'ntr=',ntr
        CALL TEMOL(Par,Pari,Rads,Mass,Kpars,Nh,Its,Dts,Ntr,Ins)
c . . . обход интерпол€ции шар-трубка при фиксировании ионосферы
        IF ( Mass(13).NE.0 ) CALL INST(dolg,Ntsl,Nl,Ntr,Ddolgt,Kdf,
     &                                 Ldor,Isp,Par,Pari,Pole,Nr,Ni,
     &                                 Park,Ks,Int_,Rads,Nh,Its,Dtets,
     &                                 Kpars,Qom,Qmax,Iqo,Mast)
        nfile = 9
        Pgl1(6:Kpars,1:Nh,1:Its,j) = Par(6:Kpars,1:Nh,1:Its)
        Parj(1:Nh,1:Its,j) = Parj1(1:Nh,1:Its)
      ENDDO ! j
!------------------------------------------------------------------------------
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of longitude
!------------------------------------------------------------------------------
      key = 1
! pole smoothing
      CALL BONPGL_NP(Pgl1,Kpars,Nh,Its,Ids,18)
      CALL BONPGL_NP(Pgl1,Kpars,Nh,Its,Ids,19)
      CALL BOTITE(Pgl1,Kpars,Nh,Its,Ids)
      END SUBROUTINE GLQINT_BAS
!------------------------------------------------------------------------------
      subroutine bonPGL_np(pgl,kpars,nh,its,ids,np)
	! smoothing in polar for np- number parametr
      dimension pgl(kpars,nh,its,ids)
      i2=its-1
      do k=1,nh
c      ssp - sum s.pole
c      snp - sum n.pole
        s np=0.
        s sp=0.
        do 1 j = 1 , ids
         snp=snp+pgl(np,k,2,j)
         ssp=ssp+pgl(np,k,i2,j)
   1    continue
c
        unp=snp/ids
        usp=ssp/ids
        do  j=1,ids
          pgl(np,k,1,j)=unp
          pgl(np,k,its,j)=usp
        end do 
    
      end do
      return
      end
