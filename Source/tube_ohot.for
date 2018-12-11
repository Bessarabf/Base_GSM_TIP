!  ver 18.04.14 - add nv to interface
!  Hot O influence
      subroutine tube_OH(ne,nx,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,vio1,vih1,
     *                vihe1,dolm,qomt,qmax,iqo,mast,mass,cOhot,tOhot,nv)
      
      dimension ht(nv),tt(nv),co2(nv),cn2(nv),co(nv),ch(nv),
     *          che(nv),tn(nv),vnq(nv),vnu(nv),vnv(nv),cim(nv),
     *          cio(nv),cih(nv),cihe(nv),vio(nv),vih(nv),
     *          vihe(nv),ti(nv),te(nv),vdu(nv),vdv(nv),cio1(nv),
     *          cih1(nv),cihe1(nv),ti1(nv),te1(nv),qo(nv),qsm(nv),
     *          vio1(nv),vih1(nv),vihe1(nv),cOhot(nv),tOhot(nv),
     *          mast(40),mass(30)

    
	real,allocatable:: lam(:),lamp(:),lamm(:) 
      allocatable al(:),alp(:),alm(:),bet(:),beta(:), 
     *            betap(:),betam(:),ga(:),gap(:),gam(:),gamm(:),
     *            delu(:), tempi(:),tempe(:),col(:)
      allocate (al(nv),alp(nv),alm(nv),bet(nv),beta(nv), 
     *          betap(nv),betam(nv),ga(nv),gap(nv),gam(nv),gamm(nv),
     *          delu(nv),tempi(nv),tempe(nv),col(nv))
 
	allocate(lam(nv),lamp(nv),lamm(nv))

      data dte/1.e-1/
!	,bkt/1.3807e-16/,aem/1.66057e-24/
c    *   ,pi/3.14159265359/
  100 format('    tube',i7)
     
	i4=nx/2
      i5=(nx+1)/2
      i1=1
      k1=i1
      k2=i4
      k3=i5+1
	IF(i4.eq.i5) THEN
	   i2=i4
	ELSE
	   i2=nx
	END IF

c      if(i4.eq.i5)i2=i4
c      if(i4.ne.i5)i2=nx
c
c     pi18=pi/180.
c

   1  continue

c
c     do i=i1,i2
c       if(ne.eq.1)cio(i)=1.e00*exp(-abs(ht(i)*1.e-5-300.)/400.)*
c    *     (1.e0+1.e0*cos(tt(i)))*(1.e3+1.e3*sin(dolm*pi18))*
c    *     (utt/3600.)**2
c       if(ne.eq.2)cih(i)=1.e00*exp(-abs(ht(i)*1.e-5-1000.)/400.)*
c    *     (1.e1+1.e1*cos(tt(i)*2.))*(1.e2+1.e2*cos(dolm*pi18))*
c    *     (utt/3600.)**2
c       if(ne.eq.3)cihe(i)=1.e-3
c       if(ne.eq.4)ti(i)=tn(i)
c       if(ne.eq.5)te(i)=tn(i)
c     end do
c     goto 25
c
      if(ne.le.3)then
	
         call alga_OH(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *                ti,te,vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,
     *                qsm,al,ga,qomt,qmax,iqo,cOhot,tOhot,mass,nv)

     	   call lambet(ne,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,vnq,cim,
     *               cio,cih,cihe,vio,vih,vihe,ti,te,col,vdu,vdv,beta,
     *               lam,vio1,vih1,vihe1,al,ga,cio1,cih1,cihe1,dolm)
	    

      end if
      if(i1.ne.1) goto 2
        if(ne.lt.4) goto 33
          if(ne.eq.4)then
            itr=mast(10)
          else
            itr=mast(12)
          end if
          tiit=ti(i1)
          teit=te(i1)
          do 32 i=1,itr
            call algat_OH(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,cim,cio,cih,cihe,
     *                    vio,vih,vihe,tiit,teit,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,alte,gate,qomt,qmax,
     *                    iqo,cOhot,tOhot,mass,nv)
            if(ne.eq.4)then
              teit1=teit
              tiit1=tiit*(1.+dte)
            else
              tiit1=tiit
              teit1=teit*(1.+dte)
            end if
		   
            call algat_OH(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,cim,cio,cih,cihe,
     *                    vio,vih,vihe,tiit1,teit1,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,altep,gatep,qomt,qmax,
     *                    iqo,cOhot,tOhot,mass,nv)
            if(ne.eq.4)then
              teit2=teit
              tiit2=tiit*(1.-dte)
            else
              tiit2=tiit
              teit2=teit*(1.-dte)
            end if
		 
            call algat_OH(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,
     *                    cim,cio,cih,cihe,vio,vih,vihe,tiit2,
     *                    teit2,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,altem,gatem,qomt,qmax,
     *                    iqo,cOhot,tOhot,mass,nv)
            if(ne.eq.4)then
              dfte=(gatep/altep-gatem/altem)/(tiit*2.*dte)
              tiit=(gate/alte-tiit*dfte)/(1.-dfte)
              if(tiit.lt.tn(i1))tiit=tn(i1)
            else
              dfte=(gatep/altep-gatem/altem)/(teit*2.*dte)
              teit=(gate/alte-teit*dfte)/(1.-dfte)
              if(teit.lt.tn(i1))teit=tn(i1)
            end if
   32     continue
         	
          call algat_OH(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                  che,tn,vnq,vnu,vnv,cim,cio,cih,cihe,
     *                  vio,vih,vihe,tiit,teit,vdu,vdv,cio1,cih1,
     *                  cihe1,ti1,te1,qo,qsm,alte,gate,qomt,qmax,
     *                  iqo,cOhot,tOhot,mass,nv)
          al(i1)=alte
          ga(i1)=gate
          if(ne.eq.4)then
            ti(i1)=tiit
          else
            te(i1)=teit
          end if
   33   continue
        bet(i1)=0.
        gamm(i1)=ga(i1)/al(i1)
        goto4
    2 continue
        if(ne.gt.3)goto3
          bet(i1)=0.
          if(ne.eq.1)gamm(i1)=1.e-15
          if(ne.eq.2)gamm(i1)=5.e-3
c         if(ne.eq.1)gamm(i1)=1.e-5
c         if(ne.eq.2)gamm(i1)=1.e-3
          if(ne.eq.3)gamm(i1)=1.e-8
ccc       gamm(i1)=1.e-3
          goto4
    3 continue
        bet(i1)=0.
cc      gamm(i1)=tn(i1)*3.
cc      if(ne.eq.4)gamm(i1)=tn(i1)
        gamm(i1)=tn(i1)*3.
        if(ne.eq.4)gamm(i1)=tn(i1)*3.
c       gamm(i1)=tn(i1)
c       if(ne.eq.4)gamm(i1)=tn(i1)
c       gamm(i1)=tn(i1)*5.
c       if(ne.eq.4)gamm(i1)=tn(i1)*5.
c       if(ne.eq.4)gamm(i1)=tn(i1)*10.
cc      if(ne.eq.5)te(i1)=tn(i1)*3.
        if(ne.eq.4)ti(i1)=tn(i1)*3.
        if(ne.eq.5)te(i1)=tn(i1)*3.
c       if(ne.eq.4)ti(i1)=tn(i1)
c       if(ne.eq.5)te(i1)=tn(i1)
c       if(ne.eq.4)ti(i1)=tn(i1)*10.
c       if(ne.eq.4)ti(i1)=tn(i1)*5.
c       if(ne.eq.5)te(i1)=tn(i1)*5.
    4 continue
      if(ne.lt.4)goto35
        if(i2.ne.nx)goto35
          if(ne.eq.4)then
            itr=mast(10)
          else
            itr=mast(12)
          end if
          tiit=ti(i2)
          teit=te(i2)
          do 34 i=1,itr
	     
            call algat_OH(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,
     *                    cim,cio,cih,cihe,vio,vih,vihe,
     *                    tiit,teit,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,alte,gate,qomt,qmax,
     *                    iqo,cOhot,tOhot,mass,nv)
            if(ne.eq.4)then
              teit1=teit
              tiit1=tiit*(1.+dte)
            else
              tiit1=tiit
              teit1=teit*(1.+dte)
            end if
             
		  call algat_OH(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,
     *                    cim,cio,cih,cihe,vio,vih,vihe,tiit1,
     *                    teit1,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,altep,gatep,
     *                    qomt,qmax,iqo,cOhot,tOhot,mass,nv)
            if(ne.eq.4)then
              teit2=teit
              tiit2=tiit*(1.-dte)
            else
              tiit2=tiit
              teit2=teit*(1.-dte)
            end if
	      call algat_OH(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,
     *                    cim,cio,cih,cihe,vio,vih,vihe,tiit2,
     *                    teit2,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,altem,gatem,
     *                    qomt,qmax,iqo,cOhot,tOhot,mass,nv)
            if(ne.eq.4)then
              dfte=(gatep/altep-gatem/altem)/(tiit*2.*dte)
              tiit=(gate/alte-tiit*dfte)/(1.-dfte)
              if(tiit.lt.tn(i2))tiit=tn(i2)
            else
              dfte=(gatep/altep-gatem/altem)/(teit*2.*dte)
              teit=(gate/alte-teit*dfte)/(1.-dfte)
              if(teit.lt.tn(i2))teit=tn(i2)
            end if
   34     continue
          call algat_OH(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                  che,tn,vnq,vnu,vnv,
     *                  cim,cio,cih,cihe,vio,vih,vihe,tiit,teit,
     *                  vdu,vdv,cio1,cih1,
     *                  cihe1,ti1,te1,qo,qsm,alte,gate,
     *                  qomt,qmax,iqo,cOhot,tOhot,mass,nv)
          
		al(i2)=alte
          ga(i2)=gate
          if(ne.eq.4)then
            ti(i2)=tiit
          else
            te(i2)=teit
          end if
   35 continue
      if(ne.le.3)
     *   call forvfv(ne,i1,i2,ht,tt,cim,cio,cih,cihe,
     *               alf,beta,bet,gamm,ga,al,lam,vio,vih,vihe,col)
      goto(5,8,11,14,16),ne
    5 continue
        if(i2.ne.nx)goto6
          cio(nx)=ga(nx)/al(nx)
          goto7
    6   continue
ccc       cio(i2)=1.e-3
          cio(i2)=1.e-15
c         cio(i2)=1.e-5
    7   continue
        lam(i2)=(gamm(i2)-alf*cio(i2))/bet(i2)
        goto18
    8 continue
        if(i2.ne.nx)goto9
          cih(nx)=ga(nx)/al(nx)
          goto10
    9   continue
ccc       cih(i2)=1.e-3
          cih(i2)=5.e-3
c         cih(i2)=1.e-3
   10   continue
        lam(i2)=(gamm(i2)-alf*cih(i2))/bet(i2)
        goto18
   11 continue
        if(i2.ne.nx)goto12
          cihe(nx)=ga(nx)/al(nx)
          goto13
   12   continue
ccc       cihe(i2)=1.e-3
          cihe(i2)=1.e-8
   13   continue
        lam(i2)=(gamm(i2)-alf*cihe(i2))/bet(i2)
        goto18
   14 continue
        if(i2.ne.nx) goto 15
          goto 18
   15   continue
c         ti(i2)=tn(i2)
          ti(i2)=tn(i2)*3.
c         ti(i2)=tn(i2)
c         ti(i2)=tn(i2)*10.
c         ti(i2)=tn(i2)*5.
          goto 18
   16 continue
        if(i2.ne.nx) goto 17
          goto18
   17   continue
cc        te(i2)=tn(i2)*3.
          te(i2)=tn(i2)*3.
c         te(i2)=tn(i2)
c         te(i2)=tn(i2)*5.
   18 continue
      if(ne.le.3)
     *  call backfv(ne,i1,i2,ht,tt,cim,cio,cih,cihe,ti,te,
     *              bet,gamm,ga,al,lam)
      if(ne.gt.3)then
        delu(i1)=0.
        delu(i2)=0.
        if(ne.eq.4)then
          itr=mast(10)
        else
          itr=mast(12)
        end if
        do40i=1,itr
          bet(i1)=0.
          gamm(i1)=delu(i1)
	
          call alga_OH(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *             vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *             ti,te,vdu,vdv,cio1,cih1,
     *             cihe1,ti1,te1,qo,qsm,al,ga,qomt,qmax,
     *             iqo,cOhot,tOhot,mass,nv)
	  
          call lambet(ne,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,vnq,cim,
     *            cio,cih,cihe,vio,vih,vihe,ti,te,col,vdu,vdv,beta,lam,
     *            vio1,vih1,vihe1,al,ga,cio1,cih1,cihe1,dolm)
       

          do j=i1,i2
            if(ne.eq.4)then
              tempi(j)=ti(j)*(1.+dte)
              tempe(j)=te(j)
            else
              tempi(j)=ti(j)
              tempe(j)=te(j)*(1.+dte)
            end if
          end do
      call alga_OH(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *             vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *             tempi,tempe,vdu,vdv,cio1,cih1,
     *             cihe1,ti1,te1,qo,qsm,alp,gap,qomt,qmax,iqo,
     *             cOhot,tOhot,mass,nv)

	     

      call lambet(ne,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,vnq,cim,
     *            cio,cih,cihe,vio,vih,vihe,tempi,tempe,col,vdu,vdv,
     *            betap,lamp,vio1,vih1,vihe1,alp,gap,
     *            cio1,cih1,cihe1,dolm)
	      
          do 43 j=i1,i2
            if(ne.eq.4)then
              tempi(j)=ti(j)*(1.-dte)
              tempe(j)=te(j)
            else
              tempi(j)=ti(j)
              tempe(j)=te(j)*(1.-dte)
            end if
   43     continue
      call alga_OH(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *             vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *             tempi,tempe,vdu,vdv,cio1,cih1,cihe1,
     *             ti1,te1,qo,qsm,alm,gam,qomt,qmax,iqo,
     *             cOhot,tOhot,mass,nv)

	call lambet(ne,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,vnq,cim,
     *     cio,cih,cihe,vio,vih,vihe,tempi,tempe,col,vdu,vdv,betam,
     *     lamm,vio1,vih1,
     *     vihe1,alm,gam,cio1,cih1,cihe1,dolm)
	      

      if(ne.eq.4)then
        call forvt(ne,i1,i2,ht,tt,cim,cio,cih,cihe,beta,
     *  bet,gamm,ga,al,lam,betap,gap,alp,lamp,betam,gam,alm,lamm,
     *  ti)
      else
        call forvt(ne,i1,i2,ht,tt,cim,cio,cih,cihe,beta,
     *  bet,gamm,ga,al,lam,betap,gap,alp,lamp,betam,gam,alm,lamm,
     *  te)
      end if
      call backt(i1,i2,delu,bet,gamm)
          do 44 j=i1,i2
            if(ne.eq.4)then
              ti(j)=ti(j)+delu(j)
            if(ti(j).lt.tn(j))ti(j)=tn(j)
            if(ti(j).gt.1.e4)ti(j)=1.e4
            else
              te(j)=te(j)+delu(j)
            if(te(j).lt.tn(j))te(j)=tn(j)
            if(te(j).gt.1.e4)te(j)=1.e4
            end if
   44     continue
   40 continue
      end if
      if(ne.gt.3)goto31
        if(i1.eq.1)lam(i1)=0.
        i3=i1+1
        if(i1.ne.1)lam(i1)=lam(i3)*bdip(ht(i1),tt(i1))
        do20m=i3,i2
          mp=m+1
          if(m.ne.i2)goto19
            if(i2.eq.nx)lam(i2)=0.
            if(i2.ne.nx)lam(i2)=lam(i2)*bdip(ht(i2),tt(i2))
            goto20
   19     continue
          w=(lam(m)+lam(mp))*.5
          bm=bdip(ht(m),tt(m))
          lam(m)=w*bm
   20   continue
  103 format(' ',5x,'tube: Print ')
      bm=1.e6
      bm1=-1.e6
   31 continue
        do24m=i1,i2
          goto(21,22,23,29,30),ne
   21     continue
            w=lam(m)/cio(m)
c           cio(m)=cio(m)*.1
            if(m.le.k2.and.w.gt.bm)w=bm
            if(m.ge.k3.and.w.lt.bm1)w=bm1
c           cs=sqrt(bkt/(16.*aem)*(ti(m)+cio(m)/(cio(m)+
c    *         cih(m)+cihe(m))*te(m)))*2.
c           if(abs(w).gt.cs)w=sign(cs,w)
            vio(m)=w
            if(cio(m).lt.0.)vio(m)=-vio(m)
cc          if(cio(m).lt.1.e-3)cio(m)=1.e-3
            if(cio(m).lt.0.)cio(m)=1.e-15
            goto24
   22     continue
            w=lam(m)/cih(m)
c           cih(m)=cih(m)*.1
            if(m.le.k2.and.w.gt.bm)w=bm
            if(m.ge.k3.and.w.lt.bm1)w=bm1
c           cs=sqrt(bkt/(1.*aem)*(ti(m)+cih(m)/(cio(m)+
c    *         cih(m)+cihe(m))*te(m)))*2.
c           if(abs(w).gt.cs)w=sign(cs,w)
            vih(m)=w
            if(cih(m).lt.0.)vih(m)=-vih(m)
cc          if(cih(m).lt.1.e-3)cih(m)=1.e-3
            if(cih(m).lt.0.)cih(m)=1.e-15
            goto24
   23     continue
            w=lam(m)/cihe(m)
            if(m.le.k2.and.w.gt.bm)w=bm
            if(m.ge.k3.and.w.lt.bm1)w=bm1
c           cs=sqrt(bkt/(4.*aem)*(ti(m)+cihe(m)/(cio(m)+
c    *         cih(m)+cihe(m))*te(m)))*2.
c           if(abs(w).gt.cs)w=sign(cs,w)
            vihe(m)=w
            if(cihe(m).lt.0.)vihe(m)=-vihe(m)
cc          if(cihe(m).lt.1.e-3)cihe(m)=1.e-3
            if(cihe(m).lt.0.)cihe(m)=1.e-15
            goto24
   29     continue
            goto24
   30     continue
   24   continue
  101   format(' ',5x,'tube:   Print ')
  102   format(' ',10g12.3)
   25 continue
      if(i2.eq.nx) goto26
        i1=i4+1
        i2=nx
        go to 1
   26 continue
      deallocate (lam,lamp,lamm)

      deallocate (al,alp,alm,bet,beta, 
     *          betap,betam,ga,gap,gam,gamm,
     *          delu,tempi,tempe,col) 
      return
      end


