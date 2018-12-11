      subroutine trubka_OH(ntsl,nl,pari,ni,park,ks,vdr,par1,nr,mast,
     *                     dolm,nv,ddolgt,kdf,ldor,isp,pole,kpart,
     *                     cio1,cih1,cihe1,vio1,vih1,vihe1,ti1,te1,
     *                     cio,cih,cihe,vio,vih,vihe,ti,te,co2,cn2,
     *                     co,ch,che,tn,qo,qsm,cim,vnq,vnu,vnv,ht,tt,
     *                     vdv,vdu,dtt,sole,solen,nse,int,
     *                     qom,qmax,iqo,utt,cOhot,tOhot,mass)
      dimension ntsl(nl),pari(ni),park(ks),vdr(ks),par1(nr),
     *          kdf(20),cio1(nv),cih1(nv),cihe1(nv),vio1(nv),vih1(nv),
     *          vihe1(nv),ti1(nv),te1(nv),cio(nv),cih(nv),cihe(nv),
     *          vio(nv),vih(nv),vihe(nv),ti(nv),te(nv),co2(nv),
     *          cn2(nv),co(nv),ch(nv),che(nv),tn(nv),qo(nv),qsm(nv),
     *          cim(nv),vnq(nv),vnu(nv),vnv(nv),ht(nv),tt(nv),
     *          vdv(nv),vdu(nv),sole(nse),solen(nse),sih(15),sihe(15),
     *          pole(ldor/4),qom(nl),cOhot(nv),tOhot(nv),
     *          mass(30),mast(40)
      logical readfl
c . . . Сечения
      data sih/0.,0.,8.1e-9,5.e-9,2.3e-9,.9e-9,.4e-9,.2e-9,.1e-9,.1e-9,
     *         0.,0.,0.,0.,0./,
     *    sihe/0.,0.,0.,0.,3.e-9,3.4e-9,1.6e-9,3.e-9,5.e-9,5.e-9,1.7e-9,
     *        .5e-9,0.,0.,0./
      q0h=0.
      q0he=0.
c      print 900,dolm
  900 format('+trubka  dolm=',f4.0)
      do i=1,nse
        q0h=q0h+sih(i)*sole(i)
        q0he=q0he+sihe(i)*sole(i)
      end do
      do 15 j=1,nl
        nx=ntsl(j)
        call select_OH(1,j,kpart,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1
     *             ,ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm,ht,tt,vdv,vdu,par1,nr,cOhot,tOhot,NV)
        call select_OH(2,j,int,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *              ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm, ht,tt,vdv,vdu,pari,ni,cOhot,tOhot,nv)
        call select_OH(3,j,2,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *              ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm, ht,tt,vdv,vdu,park,ks,cOhot,tOhot,NV)
        call select_OH(4,j,2,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,vihe1,
     *              ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *              qo,qsm, ht,tt,vdv,vdu,vdr,ks,cOhot,tOhot,NV)
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
              call tube_oh(1,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,vnu,
     *                  vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                  vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                  vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                  cOhot,tOhot,nv)

            end do
            if(mast(4).eq.0)goto6
              lh=mast(5)
              do ith=1,lh
                call tube_oh(2,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                    vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                    vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                    vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                    cOhot,tOhot,nv)
              end do
              if(mast(6).eq.0)goto6
              lhe=mast(7)
              do ithe=1,lhe
                call tube_oh(3,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                    vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                    vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,vio1,
     *                    vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *					cOhot,tOhot,nv)
              end do
              go to 7
    6       continue
            call aprhe(nx,cihe,vihe)
            if(mast(4).eq.0)call aprh(nx,cih,vih)
    7     continue
          lt=mast(8)
          do itt=1,lt
            if(mast(9).ne.0) then
               call tube_oh(4,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                    vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                    vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                    vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                    cOhot,tOhot,nv)
            else
              call aprti(nx,tn,ti)
            end if
            if(mast(11).ne.0) then
                call tube_oh(5,nx,dtt,ht,tt,co2,cn2,co,ch,che,tn,vnq,
     *                    vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                    vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                    vio1,vih1,vihe1,dolm,qomt,qmax,iqo,mast,mass,
     *                    cOhot,tOhot,nv)
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
c$$$$$$$$
c     nnn=17
c     msumm=0
c     nnn1=nnn-1
c     do 774 kkk=1,nnn1
c       msumm=msumm+ntsl(kkk)
c 774 continue
c     kkk1=msumm*kpart
c     kkk2=kkk1+ntsl(nnn)*kpart
c     kkk1=kkk1+1
c     do 909 iii=1,9
c     par1(kkk1+7)=par1(kkk1+6)
c     kkk1=kkk1+8
c 909 continue
c$$$$$$$$$$$
      md=1
      call wwt(readfl,nfile,kpar,dolm,ddolgt,nomsl,ntsl,nl,
     *         kdf,ldor,isp,md,par1,pole,nr,mast)
      return
      end
!
