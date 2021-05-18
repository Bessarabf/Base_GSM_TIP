! subroutine cyclt2_bas, fillin, select, vdrift, trubka, aprh, aprhe, aprte, 
!      aprti
c    VER _04_2014 DINAMIC MASSIVES

      subroutine cyclt2_bas(mast,ntsl,nl,par,nr,pari,ni,par1,PAR2,b,c,
     *                  its,park,ks,gins,rads,nh,ddolgt,dtett,ddolgs,
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
          call select(ne,j,m,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *                ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *                qo,qsm,ht,tt,vdv,vdu,park,ks,nv)
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

      call wwt(readfl,nfile,kpart,dolm,ddolgt,nomsl,ntsl,nl,kdf,ldor,
     *         isp,md,par1,pole,nr,mast)


    5  continue
        nfile=8
        kpar=int
        call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
     *           kdf,ldor,isp,md,pari,pole,ni,mast)
        call trubka(ntsl,nl,pari,ni,park,ks,vdr,par1,nr,mast,dolm,nv,
     *              ddolgt,kdf,ldor,isp,pole,kpart,cio1,cih1,cihe1,
     *              vio1,vih1,vihe1,ti1,te1,cio,cih,cihe,vio,vih,vihe,
     *              ti,te,co2,cn2,co,ch,che,tn,qo,qsm,cim,vnq,vnu,vnv,
     *              ht,tt,vdv,vdu,dtt,sole,solen,nse,int,qom,qmax,iqo
     *              ,utt,mass)

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
     *          vdu,qsm)
      return
      end
!------------------------------------------------------------------------------
      subroutine fillin(ne,j,m,ntsl,nl,cio,cih,cihe,vio,vih,vihe,
     *           ti,te,vdv,vdu,par,nr)
      dimension ntsl(nl),cio(*),cih(*),cihe(*),vio(*),
     *          vih(*),vihe(*),ti(*),te(*),vdv(*),vdu(*),
     *          par(nr)
      i1=0
      if(j.eq.1)goto2
        i2=j-1
        do1i=1,i2
          i1=i1+ntsl(i)
    1   continue
    2 continue
      n=ntsl(j)
      do 5 i=1,n
        k=(i1+i-1)*m+1
        goto(3,4),ne
    3   continue
          par(k)=cio(i)
          par(k+1)=cih(i)
          par(k+2)=cihe(i)
          par(k+3)=vio(i)
          par(k+4)=vih(i)
          par(k+5)=vihe(i)
          par(k+6)=ti(i)
          par(k+7)=te(i)
          goto5
    4   continue
          par(k)=vdv(i)
          par(k+1)=vdu(i)
    5 continue
      return
      end
!------------------------------------------------------------------------------
!     VER 18.04.14 add nv to interface      
	subroutine select(ne,j,m,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,
     *           vihe1,ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *           qo,qsm,ht,tt,vdv,vdu,par,nr,NV)

      dimension ntsl(nl),par(nr),cio1(nv),cih1(nv),cihe1(nv),vio1(nv),
     *          qsm(nv),
     *          vih1(nv),vihe1(nv),ti1(nv),te1(nv),co2(nv),cn2(nv),
     *          co(nv),ch(nv),che(nv),cim(nv),tn(nv),
     *          vnq(nv),vnu(nv),vnv(nv),qo(nv),ht(nv),tt(nv),
     *          vdv(nv),vdu(nv)
      data pi/3.14159265/
  900 format(' ',10g12.4)
      cr=pi/180.
      i1=0
      if(j.NE.1) THEN
        i2=j-1
        do i=1,i2
          i1=i1+ntsl(i)
        end do
      END IF
      n=ntsl(j)
      do 7 i=1,n
        k=(i1+i-1)*m+1
        goto(3,4,5,6),ne
    3   continue
          cio1(i)=par(k)
          cih1(i)=par(k+1)
          cihe1(i)=par(k+2)
          vio1(i)=par(k+3)
          vih1(i)=par(k+4)
          vihe1(i)=par(k+5)
          ti1(i)=par(k+6)
          te1(i)=par(k+7)
!	print*,n,i,ti1(i),te1(i)
          goto7
    4   continue
          co2(i)=par(k)
          cn2(i)=par(k+1)
          co(i)=par(k+2)
          cim(i)=par(k+3)
          qo(i)=par(k+4)
          ch(i)=par(k+5)
          che(i)=par(k+6)
          tn(i)=par(k+7)
          vnq(i)=par(k+8)
          vnu(i)=par(k+9)
          vnv(i)=par(k+10)
          qsm(i)=par(k+11)
          
          goto7
    5   continue
          ht(i)=par(k)
          tt(i)=par(k+1)*cr
          goto7
    6   continue
          vdv(i)=par(k)
          vdu(i)=par(k+1)
    7 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine vdrift(i,j,idt,ntsl,nl,ddolg,dtet,fi,ht,tt,
     *           potef,vdu,vdv,ntr,nl2)

      dimension tt(*),ntsl(nl),vdu(*),vdv(*),
     *          potef(ntr,idt,nl2),ht(*)

      data re/6371.02e5/,pi/3.14159265359/
      im=i-1
      if(i.eq.1)im=idt
      ip=i+1
      if(i.eq.idt)ip=1
      dt=dtet/90.*pi
      dd=ddolg/90.*pi
      k2=(ntsl(j)+1)/2
      r=re+ht(1)
      tp=tt(1)
      k1=1
      jp=nl2-j
      jm=jp+1
      jd=jp-1
      st=sin(tp)
      ct=cos(tp)
    1 continue
      dff=(potef(16,ip,jp)-potef(16,im,jp))/dd
      dfu=(potef(16,i,jm)-potef(16,i,jd))/dt
      if(jp.eq.4.or.jp.eq.33)dfu=(potef(16,i,jp)-potef(16,i,jd))/dt*2.
      if(jp.eq.5.or.jp.eq.34)dfu=(potef(16,i,jm)-potef(16,i,jp))/dt*2.
      dfu=dfu*.5*st**3/(r*ct)
      dff=dff*st*st/r
      do 2 k=k1,k2
        h=ht(k)
        t=tt(k)
        if(k.eq.ntsl(j))t=tt(1)
        b=bdip(h,t)
        st=sin(t)
        ct=cos(t)
        cev=1./(b*st**3)
        vdu(k)=dff*cev
        vdv(k)=-dfu*sqrt(1.+3.*ct*ct)*cev
    2 continue
      if(k1.ne.1)go to 3
         k1=k2+1
         k2=ntsl(j)
         tp=tt(1)
c        tp=tt(k2)
         jp=j+1
         jd=j
         jm=jp+1
         st=sin(tp)
         ct=-cos(tp)
c        ct= cos(tp)
         go to 1
    3 continue
c     print4,vdu(1),vdu(k2),vdv(1),vdv(k2)
c   4 format(' ',3x,'vu(1)=',g12.3,3x,'vu(n)=',g12.3,3x,
c    *'vv(1)=',g12.3,3x,'vv(n)=',g12.3)
      return
      end
!------------------------------------------------------------------------------
      subroutine trubka(ntsl,nl,pari,ni,park,ks,vdr,par1,nr,mast,
     *                     dolm,nv,ddolgt,kdf,ldor,isp,pole,kpart,
     *                     cio1,cih1,cihe1,vio1,vih1,vihe1,ti1,te1,
     *                     cio,cih,cihe,vio,vih,vihe,ti,te,co2,cn2,
     *                     co,ch,che,tn,qo,qsm,cim,vnq,vnu,vnv,ht,tt,
     *                     vdv,vdu,dtt,sole,solen,nse,int,
     *                     qom,qmax,iqo,utt,mass)
      dimension ntsl(nl),pari(ni),park(ks),vdr(ks),par1(nr),
     *          kdf(20),cio1(nv),cih1(nv),cihe1(nv),vio1(nv),vih1(nv),
     *          vihe1(nv),ti1(nv),te1(nv),cio(nv),cih(nv),cihe(nv),
     *          vio(nv),vih(nv),vihe(nv),ti(nv),te(nv),co2(nv),
     *          cn2(nv),co(nv),ch(nv),che(nv),tn(nv),qo(nv),qsm(nv),
     *          cim(nv),vnq(nv),vnu(nv),vnv(nv),ht(nv),tt(nv),
     *          vdv(nv),vdu(nv),sole(nse),solen(nse),sih(15),sihe(15),
     *          pole(ldor/4),qom(nl),
     *          mass(30),mast(40)
      logical readfl
c . . . Сечения
      data sih/0.,0.,8.1e-9,5.e-9,2.3e-9,.9e-9,.4e-9,.2e-9,.1e-9,.1e-9,
     *         0.,0.,0.,0.,0./,
     *    sihe/0.,0.,0.,0.,3.e-9,3.4e-9,1.6e-9,3.e-9,5.e-9,5.e-9,1.7e-9,
     *        .5e-9,0.,0.,0./
!      q0h=0.
!      q0he=0.
!      do i=1,nse
!        q0h=q0h+sih(i)*sole(i)
!        q0he=q0he+sihe(i)*sole(i)
!      end do
      do 15 j=1,nl
        nx=ntsl(j)
        call select(1,j,kpart,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1
     *             ,ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm,ht,tt,vdv,vdu,par1,nr,NV)
        call select(2,j,int,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *              ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm, ht,tt,vdv,vdu,pari,ni,nv)
        call select(3,j,2,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *              ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm, ht,tt,vdv,vdu,park,ks,NV)
        call select(4,j,2,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *              ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm, ht,tt,vdv,vdu,vdr,ks,NV)
        qomt=qom(j)
        do i=1,nx
          cio(i)=cio1(i)
          cih(i)=cih1(i)
          cihe(i)=cihe1(i)
          vio(i)=vio1(i)
          vih(i)=vih1(i)
          vihe(i)=vihe1(i)
          ti(i)=ti1(i)
          te(i)=te1(i)
          if(te(i).lt.tn(i)) te(i)=tn(i)
        end do
        lit=mast(1)
        do 14 itct=1,lit
          lc=mast(2)
          do 7 itc=1,lc
            lo=mast(3)
            do ito=1,lo
              call tube(1,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,vnu,
     *                  vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                  vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                  vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                  nv)

            end do
            if(mast(4).eq.0)goto6
              lh=mast(5)
              do ith=1,lh
                call tube(2,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                    vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                    vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                    vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                    nv)
              end do
              if(mast(6).eq.0)goto6
              lhe=mast(7)
              do ithe=1,lhe
                call tube(3,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                    vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                    vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,vio1,
     *                    vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,nv)
              end do
              go to 7
    6       continue
            call aprhe(nx,cihe,vihe)
            if(mast(4).eq.0)call aprh(nx,cih,vih)
    7     continue
          lt=mast(8)
          do itt=1,lt
            if(mast(9).ne.0) then
               call tube(4,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                   vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                   vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                   vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                   nv)
            else
              call aprti(nx,tn,ti)
            end if
            if(mast(11).ne.0) then
                call tube(5,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                    vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                    vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                    vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                    nv)
            else
              call aprte(nx,tn,te)
            end if
          end do
   14   continue
c       call ryt(nx,ht,tt,cim,cio,cih,cihe,vio,vih,vihe,ti,te,co2,
c    *           cn2,co,ch,che,qo,vnq,vnu,vnv,vdu,vdv,tn)
        call fillin(1,j,kpart,ntsl,nl,cio,cih,cihe,vio,vih,vihe,
     *              ti,te,vdv,vdu,par1,nr)
   15 continue
      readfl=.false.
      nfile=14
      kpar=kpart

      md=1
      call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
     *         kdf,ldor,isp,md,par1,pole,nr,mast)
c      print 900,dolm
  900 format('+trubka  dolm=',f4.0)

      return
      end
!------------------------------------------------------------------------------
      subroutine aprh(nx,cih,vih)
      dimension cih(*),vih(*)
      do 1 i=1,nx
        cih(i)=1.e-3
        vih(i)=0.
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine aprhe(nx,cihe,vihe)
      dimension cihe(*),vihe(*)
      do 1 i=1,nx
        cihe(i)=1.e-3
        vihe(i)=0.
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine aprte(nx,tn,te)
      dimension tn(*),te(*)
      do 1 i=1,nx
c       te(i)=tn(i)*5.
        te(i)=tn(i)
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine aprti(nx,tn,ti)
      
      dimension tn(*),ti(*)
      do 1 i=1,nx
c       ti(i)=tn(i)*5.
        ti(i)=tn(i)
    1 continue
      return
      end
!------------------------------------------------------------------------------
