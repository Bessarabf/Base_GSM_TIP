      subroutine inter1(k,m,ntet,x1,x2,xi,par,plm,int,
     *           kpars,nh,its)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     shar - trubka interpolation
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      dimension par(kpars,nh,its),pl(3),plm(int)
      integer intcon(6),ishcon(6),inttemp(5),ishtemp(5)
!  sootvetstvie nomerov parametrov shara & inerpol.parametrov. 
!  Concentration & ionization
      data ishcon/1,2,3,6,16,18/, intcon/1,2,3,4,5,13/ 
!  Temperatures & velocities
      data ishtemp/7,10,11,12,19/,intTEMP/8,9,10,11,14/
      if(k.ne.2) then
        nt1=m
        nt2=m+1
        its1=ntet
        its2=its1
      else
        its1=ntet
        its2=ntet+1
        nt1=m
        nt2=nt1
      end if 

      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
      is=1
      j=1
      i=is

      do knum=1,6  ! interpolation for concentration
         i=ishcon(knum)
         j=intcon(knum)
c ******
      if(par(i,nt1,its1).le.0.or.par(i,nt2,its2).le.0.)print 909,
     *      i,nt1,its1,i,nt2,its2,par(i,nt1,its1),par(i,nt2,its2)
  909 format('inter1!!!! i,nt1,its1,i,nt2,its2,par(i,nt1,its1),
     *       par(i,nt2,its2)',' - subr. inter1'/' ',6i4,2g12.3)
c*******
        pl(1)=alog(par(i,nt1,its1))
        pl(2)=alog(par(i,nt2,its2))
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(j)=exp(pl(3))
      end do 
      plm(7)=2.e6  ! He density
      do knum=1,5   ! interpolatin temperatures and velocities 
         i=ishtemp(knum)
         j=inttemp(knum)
         pl(1)=par(i,nt1,its1)
         pl(2)=par(i,nt2,its2)
         pl(3)=w1*pl(2)+w2*pl(1)
         plm(j)=pl(3)
      end do 
 
      step=28.9*plm(8)**(-0.25)  ! H - density
      plm(6)=10.**step
c     plm(6)=(10.**step)*1.e-1
        argum=par(13,nt1,its1)+par(14,nt1,its1)+par(15,nt1,its1)
        pl(1)=alog(argum)
        argum=par(13,nt2,its2)+par(14,nt2,its2)+par(15,nt2,its2)
        pl(2)=alog(argum)
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(12)=exp(pl(3))       ! sum ionization
      return
      end

