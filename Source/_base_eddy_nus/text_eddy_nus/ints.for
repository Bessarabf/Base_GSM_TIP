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
