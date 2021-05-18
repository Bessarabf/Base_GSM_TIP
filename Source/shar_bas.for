! subroutine shar_bas, globr, globrw, vmion, zamion, atomh, glqint_bas, 
!      bonPGL_np, molio2, gsmo, progp, botite, molion, molio3Nit_bas, 
!      ionize_AE, ionizv, temol, difte, hplm, inter1, inter1h, ionizu, 
!      ionize, ioniz_AE, inst
! function  erf, gamma, qc, qf, qp, pmte, pntrf
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
      call glqint_bas(pgl,par,pari,kpars,nh,its,ids,kdf,ldor,pole,
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
!------------------------------------------------------------------------------
c  version GSM  15.04.14  
	subroutine globr(pgl,par,kpars,ddolgs,dtets,nh,
     *                 kdf,ldor,isp,pole,nr,its,ids,mass)
      dimension pgl(kpars,nh,its,ids),par(kpars,nh,its)
     *         ,kdf(20),pole(ldor/4)
     *         ,mass(30)
      logical readfl
c 772 format (' globr - whod ',i3,2e15.5,5i7)
c     print 772,kpars,ddolgs,dtets,nh,ldor,nr,its,ids
      readfl=.true.
      nfile=5
      npgl=kpars*nh*its*ids

      call globrw (nfile,readfl,pgl,pole,kpars,ldor,nh,
     *            kdf,isp,npgl,its,ids,mass)
              
c 777 format (' globr - wihod ',3e13.4)
c     print 777,pgl(1,1,1,1),pgl(2,1,1,1),pgl(3,1,1,1)
      return
      end
!------------------------------------------------------------------------------
c  version GSM  15.04.14
       subroutine globrw(nfile,readfl,pgl1,pole,kpar,ldor,
     *                  nh,kdf,isp,npgl,its,ids,mass)
 
      dimension pgl1(npgl),kdf(20),mass(30),pole(ldor/4)
      logical readfl
c
  667 format(' globrw - begin')
  767 format(' globrw - end')
  700 format(' globrw :   ********  error  ********'/
     *   '  npg=',i8,'  >   npgl=',i8,'  !!!!!!  STOP  !')
c
      print 667

      npg=kpar*nh*its*ids
      if(npgl.lt.npg) go to 9
        nf=5
        if(nfile.le.8) nf=4
        isp=kdf(nfile)+1
        mdor=ldor/4
        ndor=npgl/mdor
        nost=npgl-ndor*mdor
        if(readfl)go to 4
          k=1
          do 2 j=1,ndor
            do 1 i=1,mdor
              pole(i)=pgl1(k)
              k=k+1
    1       continue
            write(nf,rec=isp)pole
            isp=isp+1
    2     continue
          read(nf,rec=isp)pole
          do 3 i=1,nost
            pole(i)=pgl1(k)
            k=k+1
    3     continue
          write(nf,rec=isp)pole
          go to 8
    4   continue
          k=1
          do 6 j=1,ndor
            read(nf,rec=isp) pole
            do 5 i=1,mdor
              pgl1(k)=pole(i)
              k=k+1
    5       continue
            isp=isp+1
    6     continue
          if(nost.eq.0)go to 11
            read(nf,rec=isp)pole
            do 7 i=1,nost
              pgl1(k)=pole(i)
              k=k+1
    7       continue
   11     continue
    8   continue
        go to 10
    9 continue
        print 700,npg,npgl
        stop
   10 continue
      print 767
      return
      end
!------------------------------------------------------------------------------
      subroutine vmion(pgl1,kpars,nh,its,rads,dtett,
     *           dtet,ddolgs,potef,ntr,nl2,ids,vim,vid,vir,mast,mass)
      dimension pgl1(kpars,nh,its,ids),mast(40),mass(30),rads(nh),
     *          vim(nh,its,ids),vid(nh,its,ids),vir(nh,its,ids),
     *          potef(ntr,ids,nl2)
      real nu0
      double precision ami,ae,ame,cr
      data nu0/1.e-9/,ami/30./,ae/1.6e-24/,e/1.6e-20/,b0/0.3/,
     *     pi/3.1415926/,re/6371.02e5/
      cr=180./pi
      ame=ae*ami/e
      nmd=360./ddolgs+0.1
      iteq=(its+1)/2
      itd=its-1
      do 12 j = 1 , ids
        do 11 itet = 2 , itd
c. . . двойной шаг для 10 - 5 гр. сетки
c        it=itet*2-1
c
      tet=dtet*(itet-1)/cr
      dip=1.
      if(tet.eq.0.)dip=pi/2.
      if(tet.eq.pi)dip=-pi/2.
      if(tet.eq.pi/2.)dip=0.
      if(dip.ne.1.)go to 5
        t=tan(tet)
        dip=atan(2./t)
    5 continue
      si=sin(dip)
      ci=cos(dip)
      sc=si*ci
      ct=cos(tet)
      st=sin(tet)
      bb=b0*sqrt(1+3.*ct*ct)
      be=ame/bb
      bnu=be*nu0
c
c 774 format (' vmion - 2 ')
c     print 774
      dd2=2*ddolgs/cr
      if(j.ne.1)go to 2
        def=(potef(ntr,2,itet)-potef(ntr,nmd,itet))/dd2
        go to 3
    2 continue
      if(j.ne.nmd)go to 4
         def=(potef(ntr,1,itet)-potef(ntr,nmd-1,itet))/dd2
         go to 3
    4 continue
c 775 format (' vmion - 3 ')
c     print 775
      def=(potef(ntr,j+1,itet)-potef(ntr,j-1,itet))/dd2
    3 continue
      defs=def/st/bb
      det=(potef(ntr,j,itet+1)-potef(ntr,j,itet-1))/2./dtett*cr
c . . . Граница замкнутых-разомкнутых
      if(itet.eq.4.or.itet.eq.33)det=(potef(ntr,j,itet)
     *                           -potef(ntr,j,itet-1))/dtett*cr
      if(itet.eq.5.or.itet.eq.34)det=(potef(ntr,j,itet+1)
     *                           -potef(ntr,j,itet))/dtett*cr
      detb=det/bb
      r=re+rads(ntr)
      et=-detb/r
      el=-defs/r
c
      if(itet.ne.iteq) go to 6
        tetd=tet-dtett/cr
        stet=1.-sin(tetd)**2
        er=-(potef(ntr,j,itet)-potef(ntr,j,itet-1))/(r*stet)
        er=er/bb
        go to 7
    6 continue
      er=-et*ci/si
    7 continue
      do 1 kk=1,nh
        g=bnu*((pgl1(1,kk,itet,j)+pgl1(2,kk,itet,j))/2.
     * +pgl1(3,kk,itet,j)/3.)
        g2=g*g
        alf=1.+g2
        etb=et/g
        etl=el/g
        etr=er/g
        a=pgl1(10,kk,itet,j)+etr
        b=pgl1(11,kk,itet,j)+etb
        c=pgl1(12,kk,itet,j)+etl
c       *****
        if (mast(13).eq.0) then
        a=pgl1(10,kk,itet,j)
        b=pgl1(11,kk,itet,j)
        c=pgl1(12,kk,itet,j)
        end if
c       *******
        if (mass(6).eq.0) then
        a=etr
        b=etb
        c=etl
        end if
c       *******
c 678 format (' vmion - cycl')
c     print 678
        vir(kk,itet,j)=((si*si+g2)*a+sc*b+g*ci*c)/alf
        vim(kk,itet,j)=((ci*ci+g2)*b+sc*a-g*si*c)/alf
        vid(kk,itet,j)=(g*si*b-g*ci*a+g2*c)/alf

    1 continue
   11 continue
   12 continue

c 900 format(' ',10g12.3)
c 901   format(' vir ')
c 902   format(' vim ')
c 903   format(' vid ')
c 777 format (' vmion - wihod ')
c     print 777
      return
      end
!------------------------------------------------------------------------------
      subroutine zamion(gins,ins,its,ids,nh,uts,mass)

      dimension gins(ins,nh,its,ids),mass(30)
	allocatable s(:)
	allocate (s(ids))
      if(mass(13).ne.0)go to 1
c     . . .
        iutss=nint(uts/3600.)
c       iutss=uts/3600.
c     . . .
   11   if(iutss.lt.24) go to 10
          iutss=iutss-24
          go to 11
   10   continue
        nd=iutss-mass(14)
        do 2 np=1,ins
          do 3 k=1,nh
            do 4 i=1,its
              do 5 j=1,ids
                jn=j+nd
                if(jn.gt.ids)jn=jn-ids
                if(jn.lt.1)jn=jn+ids
                s(j)=gins(np,k,i,jn)
    5         continue
              do 6 j=1,ids
                gins(np,k,i,j)=s(j)
    6         continue
    4       continue
    3     continue
    2   continue
    1 continue
      deallocate (s)

      return
      end
!------------------------------------------------------------------------------
      subroutine atomh(pgl,kpars,nh,its,ids,gkoor,rads,
     *                 fs,fa,ap,uts,day)
c     . . . расчет атомарного водорода по MSIS до 520 км
      integer day
      dimension pgl(kpars,nh,its,ids),tm(2),gkoor(2,its,ids)
     *         ,apm(7),rads(nh)
!      dimension dm(8) ! MSIS 90
      dimension dm(9)  ! MSIS 2000
      data om/7.27e-5/,pi/3.14159/
      iyd=80*1000+day
c     vys=rads(nh)/1.e5
      do l=1,7
        apm(l)=ap
      end do
c     write(*,*) fs,fa,ap
      do 88 i = 1 , its
        do 88 j = 1 , ids
           fig=gkoor(1,i,j)
           fig=90.-fig
           dgeo=gkoor(2,i,j)/180.*pi
           dol=gkoor(2,i,j)
           tault=uts+dgeo/om
           tault=tault/3600.
  111      if(tault.gt.24.) then
            tault=tault-24.
            go to 111
           end if
           tau2=uts
           if(tau2.gt.86400.) tau2=tau2-86400.
           do k=15,nh
             vys=rads(k)*1.e-5
             ! MSIS2000
             call gtd7(iyd,tau2,vys,fig,dol,tault,fs,fa,apm,
     *                  1,dm,tm)
             pgl(5,k,i,j)=dm(7)
           end do
88       continue
      close(25)
      return
      end
!------------------------------------------------------------------------------
! Base version glqint 11/03-2019
!
      subroutine glqint_bas(pgl1,par,pari,kpars,nh,its,ids,kdf,ldor,
     *                  pole,isp,ddolgs,dtets,sole,solen,rads,gkoor,
     *                  uts,delta,nse,gins,ins,mass,dts,ps,park,
     *                  nr,ni,kpart,ks,ddolgt,int,ntsl,nl,parj,parj1,
     *                  ntr,qom,qmax,iqo,mast,vir,E0,FAE)
      dimension pgl1(kpars,nh,its,ids),gins(ins,nh,its,ids)
     *         ,par(kpars,nh,its),pari(ins,nh,its)
     *         ,sole(nse),solen(nse),rads(nh),gkoor(2,its,ids)
     *         ,parj(nh,its,ids),parj1(nh,its),ps(10)
     *         ,mass(30),kdf(20),pole(ldor/4),ntsl(nl),park(ks)
     *         ,qom(nl),mast(40)
      dimension E0(its,ids),FAE(its,ids)
      dimension vir(nh,its,ids)
      logical readfl
      
      data key/0/

	readfl=.false.
      if(mass(20).eq.1) then
         call molio2(pgl1,gins,rads,vir,kpars,nh,its,dts,ids,
     *                ddolgs,ins,dtets,ntr)
      end if
  	print *,' glqint '
!!! longitude cycle
      do 1 j = 1 , ids
       dolg=ddolgs*(j-1)
       
       DO i = 1 , its
        do k = 1 , nh
         do np = 1 , kpars
           !!!form 2d massiv (for compatibility with old version GSM TIP)
           par(np,k,i)=pgl1(np,k,i,j)
         end do
         do nin = 1 ,ins
          pari(nin,k,i)=gins(nin,k,i,j)
         end do
        end do
       END DO
       call ioniz_AE(sole,solen,rads,par,parj1,gkoor,uts,dtets,dolg,
     *               ddolgs,delta,nh,its,ids,nse,kpars,mass,ps,E0,FAE)
       if(mass(20).eq.0) then
         call molion (par,kpars,nh,its,pari,ins,mass,dts,ntr)
       else if(mass(20).eq.2) then
!	   call molio3S(cO2plus,cNOplus,par,ids,its,nh,kpars,pari,ins,
!     *                mass,dts,j,ntr,key)  
         call molio3nIT_bas(par,ids,its,nh,kpars,pari,ins,
     *                  mass,dts,j,ntr,key)  
       end if
       call temol (par,pari,rads,mass,kpars,nh,its,dts,ntr,ins)
c . . . обход интерполяции шар-трубка при фиксировании ионосферы
      if(mass(13).ne.0) then
         call inst (dolg,ntsl,nl,ntr,ddolgt,kdf,ldor,isp,
     *              par,pari,pole,nr,ni,park,ks,int,rads,
     *              nh,its,dtets,kpars,qom,qmax,iqo,mast)
      end if
      nfile=9
       
      do  k = 1 , nh
       do  i = 1 , its
        do  np = 6 , kpars
         pgl1(np,k,i,j)=par(np,k,i)
        end do
        parj(k,i,j)=parj1(k,i)
        end do 
       end do
    1 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of longitude 
      key=1
      ! pole smoothing
      call bonPGL_np(pgl1,kpars,nh,its,ids,18)
      call bonPGL_np(pgl1,kpars,nh,its,ids,19)
      call botite(pgl1,kpars,nh,its,ids)
      return
      end
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
!------------------------------------------------------------------------------
      subroutine molio2 (pgl,pgi,rads,vir,kpars,nh,its,dts,
     *     ids,ddolgs,ins,dtet,ntr)
c
c     расчет ион. m+ c переносом по вертикали
c

      dimension pgl (kpars,nh,its,ids),rads(nh),vir(nh,its,ids),
     *     pgi(ins,nh,its,ids)
     
      allocatable dk(:),qs(:),r1(:),r2(:),
     *     p(:),ce(:),co(:)
	

      real nu0,l1,l2,lamo
      data nu0/1.e-9/,ami/30./,ae/1.6e-24/,el/1.6e-20/,b0/0.3/,
     *     pi/3.141592/,re/6371.02e5/,bk/1.38e-16/,alfa/4.2e-7/,
     *     l1/1.6e-11/,l2/6.e-13/,ge/980.665/
  771 format(' molion2: kpars,nh,its,ids,ddolgs,dtet,ntr :'/
     *       ' ',4i6,2g12.4,i4)

      allocate (dk(nh),qs(nh),r1(nh),r2(nh),
     *     p(nh),ce(nh),co(nh))

	cr=180./pi
      ame=ae*ami/el
      gkm=ae*ami*ge/bk
      nmd=360./ddolgs+0.1
c     nd=dolm/ddolgs+1
      itd=its-1
      do 1 j=1,ids
         do 2 i=2,itd
c
c        расчет коэф. зависящих от широты
c
         tet=dtet*(i-1)/cr
         dip=1.
         if(tet.eq.0.)dip=pi/2.
         if(tet.eq.pi)dip=-pi/2.
         if(tet.eq.pi/2.)dip=0.
         if(dip.ne.1.)go to 3
           t=tan(tet)
           dip=atan(2./t)
    3 continue
         si=sin(dip)
         si2=si*si
         ct=cos(tet)
         bb=b0*sqrt(1.+3*ct*ct)
         be=ame/bb
         bnu=be*nu0
c
c        расчет коэф. диффузии и фотохимии
c
         do 4 k=1,nh
c       **************
           pr=(pgl(1,k,i,j)+pgl(2,k,i,j))/2.+pgl(3,k,i,j)/3.
           g=bnu*pr
           g2=g*g
           rj=ae*ami*nu0*pr
           dk(k)=bk*(si2+g2)/(1.+g2)/rj
           ce(k)=pgl(6,k,i,j)
           qs(k)=pgl(13,k,i,j)+pgl(14,k,i,j)+pgl(15,k,i,j)
c          if(k.ge.ntr+2) go to 5
           if(k.gt.ntr+2) go to 5
             alt=alfa*(300./pgl(9,k,i,j))
             co(k)=pgl(16,k,i,j)/(l1*pgl(1,k,i,j)+l2*pgl(2,k,i,j))
             ce(k)=pgl(6,k,i,j)+co(k)
             qs(k)=qs(k)+pgl(16,k,i,j)
             p(k)=alt*ce(k)
             go to 6
    5      continue
           alt=alfa*(300./pgi(5,k,i,j))
           ce(k)=ce(k)+pgi(1,k,i,j)
           vrn=pgl(10,k,i,j)
           vtn=pgl(11,k,i,j)
           vln=pgl(12,k,i,j)
           vior=pgi(2,k,i,j)
           viot=pgi(3,k,i,j)
           viol=pgi(4,k,i,j)
           tn=pgl(7,k,i,j)
c
           tim=tn
           if(k.ge.16)tim=pgi(6,k,i,j)
           l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
           l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
c          te=pgl(9,k,i,j)
c          l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c          l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c          l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol)
c          l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol)
c
           qs(k)=qs(k)+(l1*pgl(1,k,i,j)+l2*pgl(2,k,i,j))*
     *     pgi(1,k,i,j)
           p(k)=alt*ce(k)
c       if(p(k).lt.ry)print776,k,alt,ce(k),pgi(5,k,i,j)
c 776 format(' molion2: k,alt,ce(k),pgi(5,k,i,j)'/' ',4g13.4)
    6      continue
    4   continue
c
c       расчет коэфф. при производных
c
        n1=nh-1
        do 7 k=2,n1
        h2=rads(k+1)-rads(k-1)
        cf1=(pgl(9,k+1,i,j)+pgl(8,k+1,i,j)-pgl(9,k-1,i,j)-
     *  pgl(8,k-1,i,j))/h2
        cf2=(pgi(1,k+1,i,j)-pgi(1,k-1,i,j))/h2
        cf2=cf2*pgi(5,k,i,j)/ce(k)
        r1(k)=dk(k)*(pgl(8,k,i,j)+pgl(6,k,i,j)*pgl(9,k,i,j)/
     *  ce(k))
        r2(k)=dk(k)*(cf1+cf2+gkm)-vir(k,i,j)
    7  continue
c
c      расчет коэф. на границах
c
        h2=rads(3)-rads(1)
        r1(1)=dk(1)*(pgl(8,1,i,j)+pgl(6,1,i,j)*pgl(9,1,i,j)/ce(1))
c       cf1=pgl(9,3,i,j)+pgl(8,3,i,j)-4*(pgl(9,2,i,j)+pgl(8,2,i,j))
c       cf1=cf1+3*(pgl(9,1,i,j)+pgl(8,1,i,j))/h2
c       cf2=(pgi(1,3,i,j)-4*pgi(1,2,i,j)+3*pgi(1,1,i,j))/h2
c       cf2=cf2*pgi(5,1,i,j)/ce(1)
c       r2(1)=dk(1)*(cf1+cf2+gkm)-vir(1,i,j)
        r2(1)=0.
c       h2=rads(nh)-rads(nh-2)
        h1=rads(nh)-rads(nh-1)
        r1(nh)=dk(nh)*(pgl(8,nh,i,j)+pgl(6,nh,i,j)*pgl(9,nh,i,j)/
     *  ce(nh))
c       cf1=(pgl(9,nh,i,j)+pgl(8,nh,i,j))*3+pgl(9,nh-2,i,j)+
c    *  pgl(8,nh-2,i,j)-4*(pgl(9,nh-1,i,j)+pgl(8,nh-1,i,j))
c       cf1=cf1/h2
c       cf2=(pgi(1,nh-2,i,j)-4*pgi(1,nh-1,i,j)+3*pgi(1,nh,i,j))/h2
c       cf2=cf2*pgi(5,nh,i,j)/ce(nh)
        cf1=(pgl(9,nh,i,j)+pgl(8,nh,i,j))
     *  -(pgl(9,nh-1,i,j)+pgl(8,nh-1,i,j))
        cf1=cf1/h1
        cf2=(pgi(1,nh,i,j)-pgi(1,nh-1,i,j))/h1
        cf2=cf2*pgi(5,nh,i,j)/ce(nh)
        r2(nh)=dk(nh)*(cf1+cf2+gkm)-vir(nh,i,j)

c
c       граничное условие для концентрации
c
        pgl(6,1,i,j)=(pgl(6,1,i,j)+qs(1)*dts)/(1.+p(1)*dts)
c
c       прогонка
c
        call gsmo (pgl,kpars,nh,its,ids,rads,i,j,dts,
     *  vir,r1,r2,qs,p)
    2   continue
    1   continue
        print *,'molion2 END'
	deallocate (dk,qs,r1,r2,
     *     p,ce,co)
        return
        end
!------------------------------------------------------------------------------
      subroutine gsmo(pgl,kpars,nh,its,ids,r,i,j,dt,
     *                vr1,aal,bbc,fp,p)
      dimension pgl(kpars,nh,its,ids),r(nh),vr1(nh,its,ids),
     *          aal(nh),bbc(nh),p(nh),pg(nh),fp(nh)
      allocatable  ap(:),bp(:),ep(:),ed(:),
     *             cp(:),dpp(:),pfu(:),ppo(:)
	allocate (ap(nh),bp(nh),ep(nh),ed(nh),	
     *             cp(nh),dpp(nh),pfu(nh),ppo(nh))
	data alfs/0.10/
c     data alfs/0.25/
      np=nh-1
      do 2 k=2,np
        dk=-p(k)
        ddx=(r(k+1)-r(k-1))/2.
        ap(k)=-(ddx*dk-ddx/dt)
        bp(k)=ddx*(fp(k)+(pgl(6,k,i,j)/dt))
 2    continue
      do 6 k=2,nh
          dx=r(k)-r(k-1)
          bbc2=(bbc(k)+bbc(k-1))/2.
          aal2=(aal(k-1)+aal(k))/2.
          ep(k)=(aal2/dx-bbc2/2.)
      if(abs(ep(k)).gt.0.2)goto 10
        ed(k)=aal2/dx+bbc2/2.
          cp(k)=0.
          dpp(k)=0.
          goto 6
 10       ed(k)=0.
          cp(k)=1./ep(k)
          dpp(k)=(aal2/dx+bbc2/2.)/ep(k)
 6        continue
c         **********
          ap(1)=(r(2)-r(1))/dt/2.
          bp(1)=ap(1)*pgl(6,1,i,j)
          b00=0.
          g00=pgl(6,1,i,j)
          bp(nh)=(r(nh)-r(nh-1))/dt/2.*(pgl(6,nh,i,j))
          ap(nh)=(r(nh)-r(nh-1))/dt/2.
          fnp1=0.
      call progp(ppo,pfu,ap,bp,ep,cp,dpp,nh,b00,g00,
     *ed,fnp1)
          do 100 k=1,nh
            pgl(6,k,i,j)=(pfu(k))
            if (pgl(6,k,i,j).le.1.e-6) pgl(6,k,i,j)=1.e-3
c                    Labtam
c           vr1(k,i,j)=ppo(k)/pgl(6,k,i,j)
c                    Labtam
c                    PC
            vr1(k,i,j)=ppo(k)/pgl(6,k,i,j)
c                    PC
 100      continue
c     сглаживание
          do 1 k=2,np
            pg(k)=(1.-2.*alfs)*pgl(6,k,i,j)+
     +      alfs*(pgl(6,k-1,i,j)+pgl(6,k+1,i,j))
   1      continue
          do 3 k=2,np
            pgl(6,k,i,j)=pg(k)
   3      continue
c
c         print *,'gsmo END'
      deallocate (ap,bp,ep,ed,	
     *             cp,dpp,pfu,ppo)
          return
          end
!------------------------------------------------------------------------------
      subroutine progp(ppo,pfu,ap,bp,ep,cp,dpp,nm,b00,
     *                g00,ed,fnp1)
      dimension ppo(nm),pfu(nm),ap(nm),cp(nm),ep(nm),
     *          dpp(nm),ed(nm),bp(nm)
      allocatable pb(:),pg(:)
	allocate (pb(nm),pg(nm))

      pb(1)=b00
      pg(1)=g00
      nmm=nm-1
      do 1 i=1,nmm
        if(abs(ep(i+1)).gt.0.3)go to 2
        ab=ed(i+1)+(1.-pb(i)*ep(i+1))*ap(i+1)
        pb(i+1)=(ep(i+1)*pb(i)-1.)/ab
        pg(i+1)=ep(i+1)*(pg(i)-pb(i)*bp(i))/ab
        go to 1
    2 continue
        ab=dpp(i+1)+(cp(i+1)-pb(i))*ap(i+1)
        pb(i+1)=(pb(i)-cp(i+1))/ab
        pg(i+1)=(pg(i)-pb(i)*bp(i))/ab
    1 continue
      ppo(nm)=(1.+pb(nm)*ap(nm))*
     *(fnp1+bp(nm))-pg(nm)*ap(nm)
      pfu(nm)=pg(nm)-pb(nm)*(bp(nm)+fnp1)
      do 3 j=2,nm
        i=nm-j+1
        ppo(i)=(1.+pb(i)*ap(i))*(ppo(i+1)+bp(i))
     *  -pg(i)*ap(i)
        pfu(i)=pg(i)-pb(i)*(ppo(i+1)+bp(i))
    3 continue
c     print*,'progp end'
      deallocate (pb,pg)

      return
      end
!------------------------------------------------------------------------------
      subroutine botite(pgl,kpars,nh,its,ids)
      dimension pgl(kpars,nh,its,ids)
      i2=its-1
      do 3 np = 6 , 9
        if(np.eq.7) go to 3
       do 2 k = 1 , nh
        snp=0.
        ssp=0.
        do 1 j = 1 , ids
         snp=snp+pgl(np,k,2,j)
         ssp=ssp+pgl(np,k,i2,j)
  1     continue
        unp=snp/ids
        usp=ssp/ids
        do 5 j = 1 , ids
          pgl(np,k,1,j)=unp
          pgl(np,k,its,j)=usp
  5     continue
  2    continue
  3   continue
      return
      end
!------------------------------------------------------------------------------
      subroutine molion(par,kpars,nh,its,pari,ins,mass,dts,ntr)
      dimension par(kpars,nh,its),mass(30),pari(ins,nh,its)
      real    l1,l2,moli,lamo
      data alfa/4.2e-7/,l1/1.6e-11/,l2/6.e-13/
      do 1 ig=1,its
        do 2 i=1,nh
          moli=par(6,i,ig)
      if(par(9,i,ig).le.0.1)print 900,par(9,i,ig),i,ig
  900 format(' par(9,i,ig)=',g12.4,2i5)
          alt=alfa*(300./par(9,i,ig))
          qs=par(13,i,ig)+par(14,i,ig)+par(15,i,ig)
          if (i.ge.ntr) go to 5
            qs=qs+par(16,i,ig)
            go to 6
    5     continue
          vrn=par(10,i,ig)
          vtn=par(11,i,ig)
          vln=par(12,i,ig)
          vior=pari(2,i,ig)
          viot=pari(3,i,ig)
          viol=pari(4,i,ig)
          tn=par(7,i,ig)
c
          tim=tn
          if(i.ge.16)tim=pari(6,i,ig)
          l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
          l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
c         te=par(9,i,ig)
c         l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c         l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c         l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol)
c         l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol)
c
          qs=qs+(l1*par(1,i,ig)+l2*par(2,i,ig))*pari(1,i,ig)
    6     continue
          m=1
    4     if(m.gt.mass(11)) go to 3
            r1=alt*par(6,i,ig)**2
            r1=(qs+r1)*dts
            r2=(2*par(6,i,ig)+pari(1,i,ig))*alt
            r2=1+r2*dts
      if(r2.eq.0)print 901,par(6,i,ig),i,ig,pari(1,i,ig)
  901 format('r2=0',10g12.4)
            par(6,i,ig)=(moli+r1)/r2
            m=m+1
            go to 4
    3     continue
    2   continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
!!!!  ver 11.03.2019 - O2+ and NO+ write in file4
!!!!  NEWTON ITERATION
      subroutine molio3Nit_bas(par,ids,its,nh,kpars,pari,ins,
     *                  mass,dts,j,ntr,key) 
	                 
      dimension par(kpars,nh,its),
     *          pari(ins,nh,its),mass(30),al(10)

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!TABLE. The list of chemical reactions.
!1.		$O_{2}^{+}+e\to O+O$
!2.		$O_2^++NO \to NO^++O_2$
!3.		$O_{2}^{+}+N\to N{{O}^{+}}+O$
!4.		$O_{2}^{+}+{{N}_{2}}\to N{{O}^{+}}+NO$
!5.		$N{{O}^{+}}+e\to N+O$
!6.		${{O}^{+}}+{{O}_{2}}\to O_{2}^{+}+O$
!7.		${{O}^{+}}+{{N}_{2}}\to N{{O}^{+}}+N$
!8.		$N_{2}^{+}+O\to N{{O}^{+}}+N$
!9.		$N_{2}^{+}+{{O}_{2}}\to O_{2}^{+}+{{N}_{2}}$
!10.		$N_{2}^{+}+e\to N+N$

! reaction rates (first approcsimation)     
      
      data al/2.25e-7,6.3E-10,1.8E-10,1.0E-15,4.5e-7,
     *        2.0E-11,1.2E-12,1.4E-10,6.0E-11,4.3e-7/	
! 
! cO2I = O2+ ions on previos TIME level
! cNOI = NO+ ions on previos TIME level
! cO2pl(nh,its,ids) = massiv O2+ and 
! cNOpl(nh,its,ids) = NO 
! cn1I = O2+
! cn2I = NO+
! cn3I = O+
! cn4I = N2+
! 
! cn1I,cn2I etc - previos iteration
! cn1I_N,... - Next iteration  
! dn1I, dn2I - delta for O2+ & NO+ concentrations 
!


	do ig=1,its-1 
	  do k=1,nh
!!!!!!!!! initialisation 
            cO2=par(1,k,ig)    ! cO2	     
            cN2=par(2,k,ig)    ! cNO 
            cO=par(3,k,ig)     ! cO 
	      cNO=par(4,k,ig)    ! cNO
            if(mass(21).eq.0) then 
	      cN=0.
	    else
              cN=par(5,k,ig)
            end if
	    te=par(9,k,ig)
    ! ionization rates
	    qO2=par(13,k,ig) 
	    qN2=par(14,k,ig)
	    qNO=par(15,k,ig)
	    qO=par(16,k,ig)  
            Q3= qO  
    ! lost 
            pL3=al(6)*cO2+al(7)*cN2
    !     O+       !!!!!!!! 
           if (k.lt.ntr) then
                  cn3I=Q3/pL3     
           else 
                  cn3I=pari(1,k,ig)
           end if
           cNe=par(6,k,ig)  ! Ne значения на i-ом временном шаге 
          ! first step before iteration O2+ & NO+
 !             key=0
           if (key.eq.0) then
              cn1I=.5*par(6,k,ig)! Q1(k)/pL1(k)
              cn2I=.5*par(6,k,ig)! Q2(k)/pL2(k)
              par(18,k,ig)=cn1I  !
              par(19,k,ig)=cn2I  !
           else
              cn1I=par(18,k,ig) !!cO2pl(k,ig,j)
              cn2I=par(19,k,ig) !!!cNOpl(k,ig,j)
	   end if
                ! lost rates
           pL4=al(8)*cO+al(9)+cO2+al(10)*cNe
                ! source 
           Q4= qN2
           cn4I=Q4/pL4
           Q1= qO2+al(6)*cn3I*cO2+al(9)*cn4I*cO2
             ! Newton iteration
           dnorm=11.
           it=0
           do while(dnorm.gt.0.001)
	!     do it=1,3 
              pL1=al(1)*cNe+al(2)*cNO+al(3)*cN+al(4)*cN2
              pL2=al(5)*cNe
              Q2=qNO+(al(2)*cNO+al(3)*cN+al(4)*cN2)*cn1I+
     *           al(7)*cn3I*cN2+al(8)*cn4I*cO
              F1=cn1i*(1+pL1*dts)-par(18,k,ig)-Q1*dts
              F2=cn2i*(1+pL2*dts)-par(19,k,ig)-Q2*dts
              ! Jacobian  
              G11=1.+pL1*dts+al(1)*cn1I*dts
              G12=0.
              G21=(al(2)*cNO+al(3)*cN+al(4)*cN2)*dts
              G22=1.+pL2*dts+al(5)*cn2I*dts
              ! delta O2+ & NO+
              dn1I=-F1/G11
              dn2I=-(F2+G21*F1/G11)/G22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !dnorm=sqrt(dn1I**2+dn2I**2)
              dnorm=abs(dn1I)/cn1I
              dnorm2=abs(dn2I)/cn2I
              if (dnorm2.gt.dnorm) dnorm=dnorm2
              cn1I=cn1I+dn1I
              cn2I=cn2I+dn2I
              cNe=cn1I+cn2I+cn3I+cn4I
              it=it+1
              if(it.ge.9) then
                print*,'WARNING iterations do not converge in ',
     *                 'molion3nit. dnorm=',dnorm  
                EXIT
              end if
           end do
!!!!   correction negative values
           if(cn1I.lt.0.) cn1I=par(18,k-1,ig)
           if(cn2I.lt.0.) cn2I=par(19,k-1,ig)   
           par(18,k,ig)=cn1I
           par(19,k,ig)=cn2I
 !          write(10,*) ig,k,j,cO2pl(k,ig,j),cNOpl(k,ig,j),dnorm,it
           par(6,k,ig)=cn1I+cn2I+cn4I !!!!cO2pl(k,ig,j)+cNOpl(k,ig,j)
         end do           
	end do
	return
	end
!------------------------------------------------------------------------------
!       add to interface ddolgs 10.07.2015  
      subroutine ionize_AE(par,parj,ut,dolm,ddolgs,dtets,del,nh,its,
     *           kpars,ps,E0,FAE)
      dimension ps(10),parj(nh,its)
      dimension E0(its,*),FAE(its,*)
      real ic,ia,ip,ifo,par(kpars,nh,its)
      data  emin/100./,emax/5.e4/,h/ 1./,ba/5./
      data   pi/3.141592/
      ic=ps(1)   ! soft electron precipitation flux
      gamc=ps(2)
      e0c=ps(3)  ! characteristic energy
      ia=ps(4)   ! auroral electrons
      gama=ps(5)
      e0a=ps(6)
      ip=ps(7)
      e0p=ps(8)
      ifo=ps(9)
      e0f=ps(10)
!      jg=dolm/15.+1
      jg=dolm/ddolgs+1
c
      j=1
      phid=dolm*pi/180.
      call magsm(ut,del,phid,phism,j)
      do 1 ig=1,its
        tet=dtets*(ig-1)

        ia=FAE(ig,jg)		   ! model Pakston
        e0a=E0(ig,jg)*1.e3
	
        qfsm=qf(ifo,phism,tet)
        if(tet.ge.90.-ba.and.tet.le.90.+ba)go to 2
c
          qcsm=qc(ic,phism,tet)
          qasm=ia
c     !!! qp in south hemisphere is different!!!
          if(tet.le.90.) then
           qpsm=qp(ip,phism,tet)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           goto 3
          end if
          qpsm=qp(ip,phism,tet)
          if(dolm.ge.180.) then
            dolms=360.-dolm
          else
            dolms=dolm
          end if
          qpsm=qpsm*exp(-(dolms-10.)/30.)**2
    3   continue

c         cusp precipitation
          call ionizv(qcsm,gamc,e0c,emin,emax,h,par,parj,kpars,
     *                nh,its,ig)

c         addition precipitation
!         if(tet.ge.90.) then
!            call ionizv(qpsm,gamc,e0p,emin,emax,h,par,parj,kpars,
!     *               nh,its,ig)
!            else
c           !!different Eo=e0c for northen hemisphere!!
!            call ionizv(qpsm,gamc,e0c,emin,emax,h,par,parj,kpars,
!     *               nh,its,ig)
!          end if
c         auroral precipitation
          if(e0a.ne.0.)call ionizv(qasm,gama,e0a,emin,emax,h,par,
     *                             parj,kpars,nh,its,ig)
    2   continue
c         ground precipitation
 !         call ionizv(qfsm,gamc,e0f,emin,emax,h,par,parj,kpars,
 !    *               nh,its,ig)
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine ionizv(i0,gam,e0,emin,emax,h,par,parj,kpars,nh,its,ig)
      dimension par(kpars,nh,its),parj(nh,its)
      real i0,c(4)
      func1(e)=1.46e-15*(e**(-5.8)+5.e9*e**(-8.8))
!!    ! func2(e,e0,r)=e**(gam-0.854)*exp(-func1(e)*r**3-e/e0)
! ionization branch: O2,N2,NO,O
!      data c/1.09,1.,0.61,1./
      data c/1.09,1.,0.,1./ ! NO does not ionizate by particles 25.04.19 
      a=1.277e-13*i0
      g=gam-0.854
      g1=1.+gam
      fg=gamma(g1)
      do 1 i=1,nh
        e1=emin
        r=1.09*par(1,i,ig)+par(2,i,ig)+0.61*par(3,i,ig)
	  em=e0*(1.+8.16e-17*e0**0.44*r**1.323)
        em=0.373*em**0.1
        em=em*r**0.3
	  em=em+g*e0

	psi=g+(5.76e-14*r)*r*r*(em**(-5.8)+1.1e10*em**(-8.8))
	
        sq=1./sqrt(psi)
        r=sqrt(psi/2.)
        er=1.+erf(r)
        fi=1.253/fg*sq*er
        fff=func1(em)
        fff=fff*r*r*r
        if(fff.gt.80.)fff=80.
        qer=em**(-0.854)*exp(-fff)
        eme0=em/e0
        if(eme0.ge.80.)eme0=80.
        sum=(eme0)**g1
        sum=fi*qer*sum
      
        sumsum=exp(-eme0)
        sum=sum*sumsum
        if(sum.lt.1.e-30)sum=0.
    3   as=a*sum
c
!       parj(i,ig)=0.
        do 4 j=1,4
          jj=j
          if(j.eq.3) jj=4
          if(j.eq.4) jj=3
          parj(i,ig)=parj(i,ig)+c(jj)*par(jj,i,ig)*as
          par(j+12,i,ig)=par(j+12,i,ig)+c(jj)*par(jj,i,ig)*as
    4   continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
      function erf(x)
      double precision y,t,u,p,a1,a2,a3,a4,a5,v,w
      data p,a1,a2,a3,a4,a5/.3275911,.254829592,-.284496736,
     *1.421413741,-1.453152027,1.061405429/
      y=x
      t=1.d0/(1.d0+p*y)
      u=dexp(-y*y)
      v=t*t
      w=a1*t+a2*v
      v=v*t
      w=w+a3*v
      v=v*t
      w=w+a4*v
      v=v*t
      w=w+a5*v
      erf=1.d0-w*u
      return
      end
!------------------------------------------------------------------------------
      function gamma(x)
      gamma=1.              
      return
      end
!------------------------------------------------------------------------------
      function qc(vic,lam,teta)
      real lam,lamg,lmn,lms
      data rad/57.2957/,vc/3.0e8/
      data fmdn/15./,fmnn/25./,dfn/10./,lmn/45./,dln/25.0/
c     data fmdn/15./,fmnn/25./,dfn/10./,lmn/45./,dln/35.0/
      data fmds/15./,fmns/25./,dfs/10./,lms/00./,dls/1.e6/
c
      lamg=lam*rad
c     !for lmn=0 only!
c     if(lamg.ge.180.) lamg=360.-lamg
      tetag=teta
c
      if(tetag.ge.90.) then
        fms=(fmds+fmns)/2.+cos(lam)*(fmds-fmns)/2.
c       r=((tetag-fms)/dfs)**2
c       r=((tetag-180.+fms)/dfs)**2
        r=((tetag-175.+fms)/dfs)**2
        r=r+((lamg-lms)/dls)**2
c       qc=vic*exp(-r)*2.
c       qc=vic*exp(-r)*0.25
        qc=vic*exp(-r)*0.5
      else
        fmn=(fmdn+fmnn)/2.+cos(lam)*(fmdn-fmnn)/2.
        r=((tetag-fmn)/dfn)**2
        if(lamg.ge.(180.+lmn).and.lmn.le.180.)
     *   lamg=lamg-360.
        r=r+((lamg-lmn)/dln)**2
        qc=vic*exp(-r)
c       if(tetag.le.10.) qc=qc+vc
       end if
c
c     qc=vic*exp(-r)

c     r2=((tetag-10.)/10.)**2
c     if(tetag.ge.90.)r2=((tetag-170.)/10.)**2
c     qc=vic*exp(-r)+vc*exp(-r2)

  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qc
      return
      end
!------------------------------------------------------------------------------
      function qf(vic,lam,teta)
      real lam,lamg,lm,lmn,lms
      data rad/57.2957/,vc/0.0e8/
      data fmdn/15./,fmnn/25./,dfn/10./,lmn/195./,dln/60./
      data fmds/0./,fmns/0./,dfs/20./,lms/120./,dls/25./
c
      lamg=lam*rad
c     if(lamg.ge.180.) lamg=360.-lamg !for lm=0 only!
      tetag=teta
c
      if(tetag.ge.30.and.tetag.le.150.) then
        r=0.
        qf=vic*exp(-r)
        goto 1
        else
        qf=0.
        goto 1
      end if
c     if(tetag.ge.90.) then
c       fms=(fmds+fmns)/2.+cos(lam)*(fmds-fmns)/2.
c   !!! precipitation in South Geomagn. Anomaly!!!
c       r=((tetag-125.)/dfs)**2
cc      r=((tetag-175.+fms)/dfs)**2
c       r=r+((lamg-lms)/dls)**2
c       qf=vic*exp(-r)
c     else
c       if(lamg.ge.(180.+lmn).and.lmn.le.180.)
c    *   lamg=lamg-360.
c       if(lamg.le.(lmn-180.).and.lmn.ge.180.)
c    *   lamg=360.+lamg
c       fmn=(fmdn+fmnn)/2.+cos(lam)*(fmdn-fmnn)/2.
c       if(lamg.le.(lmn-180.).and.lmn.ge.180.)
c    *   lamg=360.+lamg
c       r=((tetag-fmn)/dfn)**2
c       r=r+((lamg-lmn)/dln)**2
c       qf=vic*exp(-r)
c      end if
c     qf=vic*exp(-r)+vc*exp(-r2)
    1 continue
  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qc
      return
      end
!------------------------------------------------------------------------------
      function qp(vic,lam,teta)
      real lam,lamg,lm,lmn,lms
      data rad/57.2957/,vc/0.0e8/
      data fmdn/15./,fmnn/20./,dfn/10./,lmn/195./,dln/60./
c     data fmdn/15./,fmnn/25./,dfn/10./,lmn/195./,dln/60./
c     data fmdn/15./,fmnn/25./,dfn/10./,lmn/165./,dln/30./
c     data fmdn/15./,fmnn/30./,dfn/10./,lmn/165./,dln/60./
c     data fmds/0./,fmns/0./,dfs/25./,lms/150./,dls/1.e6/
c     data fmds/0./,fmns/0./,dfs/20./,lms/150./,dls/25./
      data fmds/0./,fmns/0./,dfs/20./,lms/120./,dls/25./
c
      lamg=lam*rad
c     if(lamg.ge.180.) lamg=360.-lamg !for lm=0 only!
      tetag=teta
c
      if(tetag.ge.90.) then
        fms=(fmds+fmns)/2.+cos(lam)*(fmds-fmns)/2.
c   !!! precipitation in South Geomagn. Anomaly!!!
        r=((tetag-125.)/dfs)**2
cc      r=((tetag-175.+fms)/dfs)**2
        r=r+((lamg-lms)/dls)**2
        qp=vic*exp(-r)
      else
c       if(lamg.ge.(180.+lmn).and.lmn.le.180.)
c    *   lamg=lamg-360.
c       if(lamg.le.(lmn-180.).and.lmn.ge.180.)
c    *   lamg=360.+lamg
        fmn=(fmdn+fmnn)/2.+cos(lam)*(fmdn-fmnn)/2.
        if(lamg.le.(lmn-180.).and.lmn.ge.180.)
     *   lamg=360.+lamg
        r=((tetag-fmn)/dfn)**2
        r=r+((lamg-lmn)/dln)**2
c       qp=vic*exp(-r)*0.40
c       qp=vic*exp(-r)*0.30
        qp=vic*exp(-r)*0.15
       end if
c     qp=vic*exp(-r)+vc*exp(-r2)

  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qc
      return
      end
!------------------------------------------------------------------------------
      subroutine temol(par,pari,rads,mass,kpars,nh,its,dts,ntr,ins)
      dimension par(kpars,nh,its)
     *          ,mass(30),rads(nh),pari(ins,nh,its)
      real la1,la2,lamo
      data alf/1.26e-4/,na/1/,dd/0.10/
      ch=0.
      che=0.
      cih=0.
      cihe=0.
      ntr2=ntr+2
      do 1 ig=2,its-1
        do 2 i=1,ntr2
          alt=rads(i)
          tn=par(7,i,ig)
          co2=par(1,i,ig)
          cn2=par(2,i,ig)
          co=par(3,i,ig)
          qo=par(16,i,ig)
          cim=par(6,i,ig)
          te0=par(9,i,ig)
          ti=par(8,i,ig)
          qs=par(13,i,ig)+par(14,i,ig)+par(15,i,ig)
          qs=qs+qo
          pqs=qs/cim
          pqc=alf*cim
c
          tia=tn
          if(i.ge.16)tia=pari(6,i,ig)
          la1=lamo(1,tn,0.,0.,0.,tia,0.,0.,0.,tn)
          la2=lamo(2,tn,0.,0.,0.,tia,0.,0.,0.,tn)
c         la1=lamo(1,tn,0.,0.,0.,tia,0.,0.,0.,te0)
c         la2=lamo(2,tn,0.,0.,0.,tia,0.,0.,0.,te0)
c         la1=lamo(1,tn,0.,0.,0.,tn,0.,0.,0.)
c         la2=lamo(2,tn,0.,0.,0.,tn,0.,0.,0.)
c
          cio=qo/(la1*co2+la2*cn2)
          ce=cim+cio
          qe=pgfkr(qs,ce,alt,co2,cn2,co,ch,che)
          tes=te0
        if(tes.lt.0.) print22,tes,i,ig
  22    format(' tes=',e10.3,' i= ',i4,' ig=',i4)
          m=1
    4     if(m.gt.mass(16))go to 3
            p1=pnte(alt,co2,cn2,co,ch,che,tes)
            p2=pntrf(alt,co2,cn2,co,tn,tes)
            p3=pmte(cim,tes)
            p4=petd12(co2,cn2,tn,tes,tn)
            p5=petd3(co,tn,tes)
            p6=pite(cio,cih,cihe,tes)
            f1=(p1+p2+p6)*tn+p3*ti+p4+p5+pqc
            f2=p1+p2+p3+p6+pqs
            call difte(dd,alt,co2,cn2,co,ch,che,tes,tn,ti,
     *                    cim,cio,cih,cihe,dgte)
c           tes=(te0+(f1+qe)*dts)/(1.+f2*dts)
            tes=(te0+(f1+qe-tes*dgte)*dts)/(1.+(f2-dgte)*dts)
            if(tes.le.tn) tes=tn
            m=m+1
            go to 4
    3     continue
          par(9,i,ig)=tes
    2   continue
        ntr3=ntr+3
        do 5 i=ntr3,nh
          par(9,i,ig)=pari(5,i,ig)
    5   continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine difte(dd,alt,co2,cn2,co,ch,che,tes,tn,ti,
     *           cim,cio,cih,cihe,dgte)
      dte=dd*tes
      te1=tes-dte
      te2=tes+dte
            p1=pnte(alt,co2,cn2,co,ch,che,te1)
            p2=pntrf(alt,co2,cn2,co,tn,te1)
            p3=pmte(cim,te1)
            p4=petd12(co2,cn2,tn,te1,tn)
            p5=petd3(co,tn,te1)
            p6=pite(cio,cih,cihe,te1)
      gte1=-(p1+p2+p6)*(te1-tn)-p3*(te1-ti)+p4+p5
            p1=pnte(alt,co2,cn2,co,ch,che,te2)
            p2=pntrf(alt,co2,cn2,co,tn,te2)
            p3=pmte(cim,te2)
            p4=petd12(co2,cn2,tn,te2,tn)
            p5=petd3(co,tn,te2)
            p6=pite(cio,cih,cihe,te2)
      gte2=-(p1+p2+p6)*(te2-tn)-p3*(te2-ti)+p4+p5
      dgte=(gte2-gte1)/(2.*dte)
      return
      end
!------------------------------------------------------------------------------
      function pmte(cim,te)
      data c/1.91e-3/,s/-1.5/
      pmte=c*cim*te**s
      return
      end
!------------------------------------------------------------------------------
      function pntrf(alt,co2,cn2,co,tn,te)
      data c1/5.34e-10/,c2/2.71e-10/,c3/2.63e-8/,c4/7.e-5/
      pntrf1=0.
      pntrf2=0.
      if(alt.gt.1.e8)goto1
        ts=1./sqrt(te)
        pntrf1=c1*co2*ts
        pntrf2=c2*cn2*ts
    1 continue
      pntrf3=c3*co*(1.-c4*te)/tn
      pntrf=pntrf1+pntrf2+pntrf3
      return
      end
!------------------------------------------------------------------------------
      subroutine hplm(plm,hms,h,int)
      dimension plm(int),am(7)
      data am   /32.,28.,16.,30.,16.,1.,4./,re/6371.02e5/
      data alf/.8/
c     data alf/.7/
c     data alf/.5/
cc    data alf/.1/
c     data alf/.2/
c     data alf/0./
c     data alf/.95/
      bolc=1.38041e-16
      at1=1.66e-24
      g0=980.665
      ge=g0/(1+h/re)**2
      tkg=bolc*plm(8)/ge/at1
      rh=h-hms
      hr=(re+h)/(re+hms)
      i=1
    1 if(i.gt.7) go to 2
        hi=tkg/am(i)
        argum=rh/hi*hr
        if(argum.ge.100.)argum=100.
        if(i.eq.6)then
          ho=tkg/am(3)
          ho=rh/ho*hr
          plm(i)=plm(i)*((1.-alf)*exp(-argum)+alf*exp(-ho))
        else
          plm(i)=plm(i)*exp(-argum)
        end if
        i=i+1
        go to 1
    2 continue
        hi=tkg/am(4)
        argum=rh/hi*hr
        if(argum.ge.100.)argum=100.
        plm(12)=plm(12)*exp(-argum)
      ! hot O 
      tko=bolc*plm(14)/ge/at1
      hoh=tko/am(3)
      arg=rh/hoh*hr
      if(arg.ge.100.)arg=100.
      plm(13)=plm(13)*exp(-arg)
      return
      end
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
c . . . интерпол€ци€ дл€ случа€ атом. Ќ по MSIS
      subroutine inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *        kpars,nh,its)
      dimension par(kpars,nh,its),pl(3),plm(int)
      if(k.eq.2) go to 1
        nt1=m
        nt2=m+1
        its1=ntet
        its2=its1
        go to 2
    1   continue
        its1=ntet
        its2=ntet+1
        nt1=m
        nt2=nt1
    2 continue
      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
      is=1
      j=1
      i=is
    9 if(i.gt.16)go to 3
c ******
      if(par(i,nt1,its1).le.0.or.par(i,nt2,its2).le.0.)print 909,
     *      i,nt1,its1,i,nt2,its2,par(i,nt1,its1),par(i,nt2,its2)
  909 format(' i,nt1,its1,i,nt2,its2,par(i,nt1,its1),par(i,nt2,its2)',
     *     ' - subr. inter1'/' ',6i4,2g12.3)
c*******
        pl(1)=alog(par(i,nt1,its1))
        pl(2)=alog(par(i,nt2,its2))
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(j)=exp(pl(3))
        if(i.eq.3) is=3
        if(i.eq.6) is=10
        i=i+is
        j=j+1
        go to 9
    3 continue
      plm(7)=2.e6
      is=3
      i=7
      j=8
   10 if(i.gt.12)go to 6
        pl(1)=par(i,nt1,its1)
        pl(2)=par(i,nt2,its2)
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(j)=pl(3)
        if(i.eq.10) is=1
        i=i+is
        j=j+1
        go to 10
    6 continue
c     step=28.9*plm(8)**(-0.25)
c     plm(6)=10.**step
c     plm(6)=(10.**step)*1.e-1
        pl(1)=alog(par(5,nt1,its1))
        pl(2)=alog(par(5,nt2,its2))
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(6)=exp(pl(3))
c
        argum=par(13,nt1,its1)+par(14,nt1,its1)+par(15,nt1,its1)
        pl(1)=alog(argum)
        argum=par(13,nt2,its2)+par(14,nt2,its2)+par(15,nt2,its2)
        pl(2)=alog(argum)
        pl(3)=w1*pl(2)+w2*pl(1)
        plm(12)=exp(pl(3))
      return
      end
!------------------------------------------------------------------------------
      subroutine ionizu(sole,solen,rads,par,gkoor,ut,dolm,
     *       ddolgs,del,nh,its,ids,nse,kpars)
      dimension am(6),rads(nh),ai(4),sa(4,15),si(4,15),
     *       sole(nse),solen(nse),gkoor(2,its,ids),
     *       par(kpars,nh,its)
c     double precision tau0,xp,xxp,sx,reh,chep,derf,hsm,
c    *       re,pi,xlg,cr,xig,xi,cx,chab,ab,ex1,ex2
      data am/32.,28.,16.,30.,4.,1./,bk/1.38e-16/,
     *       ae/1.66e-24/,g0/980.665/,pi/3.1415926d0/,
     *       re/6371.02d5/,om/7.2722e-5/
      data si/4.63, 4.98, 2.49, 0.0,  7.0, 4.23, 6.0, 0.0 ,12.68, 9.4,
     *        7.20, 0.0 ,13.0 , 9.0,  7.5, 0.0 ,12.5, 9.36, 7.58, 0.0,
     *        13.0, 8.56, 9.12, 0.0, 15.5, 8.5 , 9.2, 0.0 ,19.53,18.00,
     *        10.5, 0.0 ,24.12,23.47,12.69,0.0 ,16.27,20.58,
     *        8.16, 0.0 , 5.93, 0.0 , 4.0 ,0.0 , 4.8,3*0.0 ,2.5,0.0,
     *        2*0.0,1.,0.0,5*0.,2.02/
      data sa/0.46, 0.34, 0.16, 0. , 1.67, 0.8, 1.00, 0. , 4.68 ,3.48,
     *        4.24, 0.  , 6.5 , 5.0, 5.0 , 0. , 8.5 , 5.4, 6.50, 0.0 ,
     *        12.0, 7.20, 8.00, 0.0, 15.5, 8.5, 9.2 , 0. ,19.00,18.50,
     *        10.5, 0.  ,24. , 24.4, 12.7, 0. ,27.3 ,27.3, 8.50, 0.0 ,
     *        16.2, 9.0 ,  4.0, 0. , 7.7 , 0.02, 0. , 0.0, 4.00, 0.0 ,
     *         0. , 0.  ,  1.6, 0.0, 0.0 , 0.0 , 0.01,0.0, 0.0 , 2.42/

      cr=180./pi
      bkg=bk/g0/ae
      k=4
      kl=nse
      nd=dolm/ddolgs+1
c
      do 1 ig=1,its
        rlat=gkoor(1,ig,nd)/cr
        rlat=pi/2.-rlat
        d olgg=gkoor(2,ig,nd)
        hl=ut+dolgg*3600./15.
        cx=sin(del) *sin(rlat)+cos(del)*cos(rlat)*cos(om*(hl-43200.))
        sx=sqrt(1.0-cx*cx)
        xi=acos(cx)
c       if(cx.lt.0.) xi=pi-xi
        xig=xi*cr
c ***** for Layman night radiation only
        if (xig.gt.90.0) then
         xlg=-0.96*cos(1.2*(xig-180.0)/cr)
        else
         xlg=0.0
        end if
c ***************************************
        do 9 i=1,k
          ai(i)=0.
    9   continue
        do 2 iv=1,nh
          i=nh-iv+1
        cone=0.
        sm=0.
          do 3 j=1,k
            cone=cone+par(j,i,ig)
            sm=sm+am(j)*par(j,i,ig)
    3     continue
          sm=sm/cone
          hsm=sm/(bkg*par(7,i,ig))
          reh=(rads(i)+re)*hsm
c ********* Chepmen function**********
          if(xig.LT.80.0) THEN
        	  chep=1.0/cx
          ELSE
	          chep=chept(reh,xi)
          END IF
          if(i.eq.nh) go to 7
            do 8 m=1,k
              ai(m)=ai(m)+(par(m,i,ig)+par(m,i+1,ig))*
     *        (rads(i+1)-rads(i))/2.
!	      if(ai(m).lt.0) print*, ai(m),m,par(m,i,ig),i,ig
    8       continue
    7     continue
          do 4 j=1,k
            fs=0.
            do 5 l=1,kl
              ab=0.
              do 6 m=1,k
                ab=ab+sa(m,l)*ai(m)
	
    6         continue
              ab=ab*1.e-18
              chab=chep*ab
c             if(chab.ge.100.0)chab=100.0
              ex1=exp(-chab)
c             if(ab.ge.100.0)ab=100.0
              ex2=exp(-ab+xlg)
              if(ex1.le.1.e-30) ex1=1.e-30
              if(ex2.le.1.e-30) ex2=1.e-30
c
              fs=fs+si(j,l)*(sole(l)*ex1+solen(l)*ex2)
    5       continue
            fs=fs*1.e-9
            mm=j+12
            if(j.eq.3) mm=16
            if(j.eq.4) mm=15
            par(mm,i,ig)=par( j,i,ig)*fs
            if(par(mm,i,ig).le.1.e-20)  par(mm,i,ig)=1.e-20
    4     continue
    2   continue
               do 11 ii = 13 , 16
              a1=alog(par(ii,nh-1,ig))
              a2=alog(par(ii,nh-2,ig))
              a3=alog(par(ii,nh-3,ig))
              ac=a3+3.*a1-3.*a2
              par(ii,nh,ig)=exp(ac)
   11      continue
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine ionize(par,parj,ut,dolm,dtets,del,nh,its,kpars,ps)
      dimension ps(10),parj(nh,its)
      real ic,ia,ip,ifo,par(kpars,nh,its)
      data  emin/100./,emax/5.e4/,h/ 1./,ba/5./
      data   pi/3.141592/
      ic=ps(1)
      gamc=ps(2)
      e0c=ps(3)
      ia=ps(4)
      gama=ps(5)
      e0a=ps(6)
      ip=ps(7)
      e0p=ps(8)
      ifo=ps(9)
      e0f=ps(10)
c
      j=1
      phid=dolm*pi/180.
      call magsm(ut,del,phid,phism,j)
      do 1 ig=1,its
        tet=dtets*(ig-1)
          qfsm=qf(ifo,phism,tet)
        if(tet.ge.90.-ba.and.tet.le.90.+ba)go to 2
c
          qcsm=qc(ic,phism,tet)
          qasm=qa(ia,phism,tet)
c     !!! qp in south hemisphere is different!!!
          if(tet.le.90.) then
c          qpsm=0.
           qpsm=qp(ip,phism,tet)
c          if(ut.ge.0..and.ut.le.25200.)then
c            qpsm=qpsm*.25*(5.+3.*sin(pi/46800.*ut-pi/26.))
c            qpsm=qpsm*2.
c          else
c            if(ut.gt.25200..and.ut.le.50400.)then
c              qpsm=qpsm*.25*(5.+3.*sin(pi/25200.*ut-pi*.5))
c            else
c              if(ut.ge.64800..and.ut.le.86400.)then
c              if(ut.ge.75600..and.ut.le.86400.)then
c                qpsm=qpsm*.25*(5.+3.*sin(pi/46800.*ut-49.*pi/26.))
c                qpsm=qpsm*3.
c              else
c                qpsm=qpsm*.5
c              end if
c            end if
c          end if
           goto 3
          end if
c         if(dolm.ge.315.or.dolm.le.60.)
c    *    qpsm=qp(ip,phism,tet)
           qpsm=qp(ip,phism,tet)
          if(dolm.ge.180.) then
            dolms=360.-dolm
          else
            dolms=dolm
          end if
          qpsm=qpsm*exp(-(dolms-10.)/30.)**2
c         qpsm=0.
    3   continue
c
c         ground precipitation
c         call ionizv(qfsm,gamc,e0f,emin,emax,h,par,parj,kpars,
c    *               nh,its,ig)
c         cusp precipitation
          call ionizv(qcsm,gamc,e0c,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)
c         addition precipitation
c         if(tet.ge.90.) then
c           call ionizv(qpsm,gamc,e0p,emin,emax,h,par,parj,kpars,
c    *               nh,its,ig)
c           else
c           !!different Eo=e0c for northen hemisphere!!
c           call ionizv(qpsm,gamc,e0c,emin,emax,h,par,parj,kpars,
c    *               nh,its,ig)
c         end if
c         auroral precipitation
          call ionizv(qasm,gama,e0a,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)
    2   continue
c         ground precipitation
          call ionizv(qfsm,gamc,e0f,emin,emax,h,par,parj,kpars,
     *               nh,its,ig)
    1 continue
c 667 format('667 - ionize')
c     print 667
      return
      end
!------------------------------------------------------------------------------
      subroutine ioniz_AE(sole,solen,rads,par,parj,gkoor,uts,dtets,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars,mass,ps,E0,FAE)
      dimension sole(nse),solen(nse),par(kpars,nh,its),parj(nh,its),
     *       mass(30),rads(nh),gkoor(2,its,ids),ps(10)
      dimension E0(its,*),FAE(its,*)
      call ionizu(sole,solen,rads,par,gkoor,uts,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars)
      if(mass(12).eq.1) then
	!  our model
         call ionize(par,parj,uts,dolm,dtets,delta,nh,its,kpars,ps)
       else if(mass(12).ge.2.and.mass(12).le.5) then
      !  our model+Zhang&Paxton or Vorobjov&Yagodkina auroral precipitation
         call ionize_AE(par,parj,uts,dolm,ddolgs,dtets,
     *                   delta,nh,its,kpars,ps,E0,FAE)
	else
	   print*,' without electron precipitation!'

      end if
      return
      end
!------------------------------------------------------------------------------
      function qa(vic,lam,teta)
      real lam,lamg,lm
      data rad/57.2957/,vc/0.0e8/
c     data fmd/15./,fmn/25./,df/5./,lm/210./,dl/40./
      data fmd/15./,fmn/25./,df/5./,lm/210./,dl/60./
c     data fmd/20./,fmn/20./,df/5./,lm/180./,dl/60./
c     data fmd/20./,fmn/20./,df/5./,lm/210./,dl/75./
c
      lamg=lam*rad
      tetag=teta
c
      fm=(fmd+fmn)/2.+cos(lam)*(fmd-fmn)/2.
      r=((tetag-fm)/df)**2
c      if(tetag.ge.90.)  r=((tetag-175.+fm)/df)**2
       if(tetag.ge.90.)  r=((tetag-180.+fm)/df)**2
c       if(lamg.ge.(180.+lm).and.lm.le.180.)
c    *   lamg=lamg-360.
        if(lamg.le.(lm-180.).and.lm.ge.180.)
     *   lamg=360.+lamg
      r=r+((lamg-lm)/dl)**2
      qa=vic*exp(-r)

c     r2=((tetag-10.)/10.)**2
c     if(tetag.ge.90.)r2=((tetag-170.)/10.)**2
c     qa=vic*exp(-r)+vc*exp(-r2)

  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qa
      return
      end
!------------------------------------------------------------------------------
      subroutine inst(dolm,ntsl,nl,ntr,ddolgt,kdf,ldor,
     *       isp,par,pari,pole,nr,ni,park,ks,int,
     *       rads,nh,its,dtets,kpars,qom,qmax,iqo,mast)
      dimension ntsl(nl),par(kpars,nh,its), pari(ni),park(ks),
     *       rads(nh),pole(ldor/4),msum(45),plm(14),qom(nl)
     *       ,kdf(20),mast(40)
      logical readfl
  900 format(' ',10g12.4)
      msum (1)=0
      do 1 nomsl=2,nl
        msum(nomsl)=msum(nomsl-1)+ntsl(nomsl-1)
    1 continue
      l=1
      do 2 nomsl=1,nl
        nsum=msum(nomsl)*int
        nt=ntsl(nomsl)
        do 3 i=1,nt
         h=park(l)
         t=park(l+1)
         tr=abs(t-90.)
c  . . .  интерполяция ниже 520 км
         IF(h.le.rads(nh)) THEN
c  . . .    интерполяция на экваторе
            if(tr.lt.0.01) then
              call find(nh,h,rads,m)
              k=1
              ntet=90./dtets+1
              xi=h
              x1=rads(m)
              x2=rads(m+1)
              if(mast(32).eq.0) then
               call inter1(k,m,ntet,x1,x2,xi,par,plm,int,
     *                     kpars,nh,its)
              else
               call inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *                      kpars,nh,its)
              end if
            else
c  . . .     интерполяция ВНЕ экватора
              m=i+ntr-1
              if(t.lt.90.) m=nt-i+ntr
              ntet=t/dtets
              xi=t
              x1=ntet*dtets
              x2=x1+dtets
              ntet=ntet+1
              k=2
              if(mast(32).eq.0) then
	        
               call inter1(k,m,ntet,x1,x2,xi,par,plm,int,
     *                     kpars,nh,its)
              else
               call inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *                      kpars,nh,its)
              end if
            end if
          ELSE
c  . . .   интерполяция выше 520 км
           m=nh
           if(tr.lt.0.01) then
c  . . .    интерполяция на экваторе
             it=its/2+1
             do lm=1,3
               plm(lm)=par(lm,nh,it)
             end do
             plm(4)=par(6,nh,it)
             plm(5)=par(16,nh,it)
             plm(7)=2.e6
             plm(8)=par(7,nh,it)
             plm(14)=par(19,nh,it)
             if(mast(32).eq.0) then
               step=28.9*plm(8)**(-0.25)
               plm(6)=10.**step
c              plm(6)=(10.**step)*1.e-1
             else
               plm(6)=par(5,nh,it)
             end if
             plm(9)=par(10,nh,it)
             plm(10)=par(11,nh,it)
             plm(11)=par(12,nh,it)
             plm(12)=par(13,nh,it)+par(14,nh,it)+par(15,nh,it)
             call hplm (plm,rads(nh),h,int)
           else
c  . . .   интерполяция ВНЕ экватора
             ntet=t/dtets
             xi=t
             x1=ntet*dtets
             x2=x1+dtets
             ntet=ntet+1
             k=2
c  . . .    интерполяция, если Н по MSIS
             if(mast(32).eq.0) then
               call inter1(k,m,ntet,x1,x2,xi,par,plm,int,
     *                     kpars,nh,its)
             else
               call inter1h(k,m,ntet,x1,x2,xi,par,plm,int,
     *                      kpars,nh,its)
             end if
             call hplm(plm,rads(nh),h,int)
           end if
         END IF
c    Сдвиг фазы меридионального ветра в Millstone Hill
cc         if(dolm.eq.0..or.dolm.eq.345.)then
cc           if(nomsl.eq.7.or.nomsl.eq.8)then
cc             plm(10)=plm(10)+6.e3
cc           end if
cc         end if
c
         call vplm (plm(9),plm(10),t)
         ll=nsum+int*(i-1)+1
         do 10 j=1,int
           pari(ll)=plm(j)
           ll=ll+1
   10    continue
         l=l+2
    3   continue
    2 continue
      j1=nl-1
      do13j=4,j1
        i1=0
        i2=j-1
        do12i=1,i2
          i1=i1+ntsl(i)
   12   continue
        i=iqo
        k=(i1+i-1)*int+1
        qo=pari(k+4)
        if(qom(j).lt.qo)qom(j)=qo
        if(qmax.lt.qo)qmax=qo
        n=ntsl(j)
        i=n-iqo+1
        k=(i1+i-1)*int+1
        qo=pari(k+4)
        if(qom(j).lt.qo)qom(j)=qo
        if(qmax.lt.qo)qmax=qo
   13 continue
      readfl=.false.
      nfile=8
      kpar=int
      md=1
      call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,
     * ntsl,nl,kdf,ldor,isp,md,pari,pole,ni,mast)
      return
      end