!  subroutine tube, alga, altv, ntnb, betnv, algat, backfv, backt, 
!      forvfv, forvt, lambet, pqji
!  function bdip, ds, dvn, bdip, lamh, lamo, pite, pitn, petd12, petd3,
!      pgfkr, piqj, nu, pnte, pntrf1, pntrf2, gst, omr, omt, la, oml, rik, sik,
!      be_lambet
!  ver 18.04.14 - add nv to interface

      subroutine tube(ne,nx,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,ti,te,
     *                vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,vio1,vih1,
     *                vihe1,dolm,qomt,qmax,iqo,mast,mass,nv)
      
      dimension ht(nv),tt(nv),co2(nv),cn2(nv),co(nv),ch(nv),
     *          che(nv),tn(nv),vnq(nv),vnu(nv),vnv(nv),cim(nv),
     *          cio(nv),cih(nv),cihe(nv),vio(nv),vih(nv),
     *          vihe(nv),ti(nv),te(nv),vdu(nv),vdv(nv),cio1(nv),
     *          cih1(nv),cihe1(nv),ti1(nv),te1(nv),qo(nv),qsm(nv),
     *          vio1(nv),vih1(nv),vihe1(nv),
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
	
         call alga(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *                ti,te,vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,
     *                qsm,al,ga,qomt,qmax,iqo,mass,nv)

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
            call algat(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                 che,tn,vnq,vnu,vnv,cim,cio,cih,cihe,
     *                 vio,vih,vihe,tiit,teit,vdu,vdv,cio1,cih1,
     *                 cihe1,ti1,te1,qo,qsm,alte,gate,qomt,qmax,
     *                 iqo,mass,nv)
            if(ne.eq.4)then
              teit1=teit
              tiit1=tiit*(1.+dte)
            else
              tiit1=tiit
              teit1=teit*(1.+dte)
            end if
		   
            call algat(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                 che,tn,vnq,vnu,vnv,cim,cio,cih,cihe,
     *                 vio,vih,vihe,tiit1,teit1,vdu,vdv,cio1,cih1,
     *                 cihe1,ti1,te1,qo,qsm,altep,gatep,qomt,qmax,
     *                 iqo,mass,nv)
            if(ne.eq.4)then
              teit2=teit
              tiit2=tiit*(1.-dte)
            else
              tiit2=tiit
              teit2=teit*(1.-dte)
            end if
		 
            call algat(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,
     *                    cim,cio,cih,cihe,vio,vih,vihe,tiit2,
     *                    teit2,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,altem,gatem,qomt,qmax,
     *                    iqo,mass,nv)
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
         	
          call algat(ne,nx,i1,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *               che,tn,vnq,vnu,vnv,cim,cio,cih,cihe,
     *               vio,vih,vihe,tiit,teit,vdu,vdv,cio1,cih1,
     *               cihe1,ti1,te1,qo,qsm,alte,gate,qomt,qmax,
     *               iqo,mass,nv)
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
	     
            call algat(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                 che,tn,vnq,vnu,vnv,
     *                 cim,cio,cih,cihe,vio,vih,vihe,
     *                 tiit,teit,vdu,vdv,cio1,cih1,
     *                 cihe1,ti1,te1,qo,qsm,alte,gate,qomt,qmax,
     *                 iqo,mass,nv)
            if(ne.eq.4)then
              teit1=teit
              tiit1=tiit*(1.+dte)
            else
              tiit1=tiit
              teit1=teit*(1.+dte)
            end if
             
		  call algat(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,
     *                    cim,cio,cih,cihe,vio,vih,vihe,tiit1,
     *                    teit1,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,altep,gatep,
     *                    qomt,qmax,iqo,mass,nv)
            if(ne.eq.4)then
              teit2=teit
              tiit2=tiit*(1.-dte)
            else
              tiit2=tiit
              teit2=teit*(1.-dte)
            end if
	      call algat(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                    che,tn,vnq,vnu,vnv,
     *                    cim,cio,cih,cihe,vio,vih,vihe,tiit2,
     *                    teit2,vdu,vdv,cio1,cih1,
     *                    cihe1,ti1,te1,qo,qsm,altem,gatem,
     *                    qomt,qmax,iqo,mass,nv)
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
          call algat(ne,nx,i2,i1,i2,dt,ht,tt,co2,cn2,co,ch,
     *                  che,tn,vnq,vnu,vnv,
     *                  cim,cio,cih,cihe,vio,vih,vihe,tiit,teit,
     *                  vdu,vdv,cio1,cih1,
     *                  cihe1,ti1,te1,qo,qsm,alte,gate,
     *                  qomt,qmax,iqo,mass,nv)
          
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
        do 40 i=1,itr
          bet(i1)=0.
          gamm(i1)=delu(i1)
	
          call alga(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *             vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *             ti,te,vdu,vdv,cio1,cih1,
     *             cihe1,ti1,te1,qo,qsm,al,ga,qomt,qmax,
     *             iqo,mass,nv)
	  
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
      call alga(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *             vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *             tempi,tempe,vdu,vdv,cio1,cih1,
     *             cihe1,ti1,te1,qo,qsm,alp,gap,qomt,qmax,iqo,
     *             mass,nv)

	     

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
      call alga(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *             vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *             tempi,tempe,vdu,vdv,cio1,cih1,cihe1,
     *             ti1,te1,qo,qsm,alm,gam,qomt,qmax,iqo,
     *             mass,nv)

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
!------------------------------------------------------------------------------
!     ver 18.04.14 add nv to interface
      subroutine alga(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *                ti,te,vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *		      al,ga,qomt,qmax,iqo,mass,NV)

      dimension ht(NV),tt(NV),co2(NV),cn2(NV),co(NV),ch(NV),
     *          che(NV),tn(NV),vnq(NV),vnu(NV),vnv(NV),cim(NV),
     *          cio(NV),cih(NV),cihe(NV),vio(NV),vih(NV),vihe(NV),
     *          ti(NV),te(NV),vdu(NV),vdv(NV),
     *          cio1(NV),cih1(NV),cihe1(NV),ti1(NV),te1(NV),
     *          qo(NV),al(NV),ga(NV),qsm(NV)
     *          ,mass(30)
      dimension hs(23),hn(23),yos(23),yon(23),tns(23),tnn(23),tes(23),
     *          tenor(23),cnes(23),cnen(23),alfa(23),tvn(23),tvs(23)
      real la1,la2,la3,la4,n2,lamh,lamo,mi,la5,la6
      data la5/1.e-9/,la6/1.4e-9/,e0/10./,qt0/2.0e14/

      if(nx.ge.45) then
        do m=1,23        !!!!!!!!!!!!!!!!!!!!!!!!!!
          m1=nx-m+1
          hs(m)=ht(m)*1.e-5
!
          hn(m)=ht(m1)*1.e-5
          yos(m)=co(m)
          yon(m)=co(m1)
          tns(m)=tn(m)
          tnn(m)=tn(m1)
          tes(m)=te(m)
          tenor(m)=te(m1)
          cnes(m)=cim(m)+cio(m)+cih(m)+cihe(m)
          cnen(m)=cim(m1)+cio(m1)+cih(m1)+cihe(m1)
        end do
        call altv(hs,yos,tns,tes,cnes,alfa,tvs,23)
        call altv(hn,yon,tnn,tenor,cnen,alfa,tvn,23)
c       print*,' nx=',nx,' i1=',i1
c       print*,' Tvs'
c       print*,tvs
c       print*,' Tvn'
c       print*,tvn
      end if
      if(ne.ne.5) goto 17
        qt=0.
        if(nx.eq.NV) goto 17
          i3=(nx+1)/2
          if(i3.le.18) goto 17
            call ntnb(nx,ht,tt,cio,cih,cihe,cim,b1,snt,snb)
            edd=alog10(snb)-11.
            ed=1.2
            if(edd.gt.0.)ed=1.3*edd**2.2+ed
c           ed=sqrt(4./3.*2.6e-12*snb)
            a=ed/e0
            qt=qt0*(1.-(1.+a)*exp(-a))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           qt=qt*1.5
ccc         qt=qt*2.0
c           qt=qt*3.0
c           if(nx.ge.63)then
c           if(nx.ge.79)then
c           if(nx.ge.93)then
c           if(nx.ge.109)then
c             qt=qt*(7.*nx-745.)/18.
c             qt=qt*(nx-107.)/2.
c линейный рост в 3 раза от 60 град широты до 70 град широты
c             qt=qt*(nx-100.)/9.
c
c линейный рост в 5 раз от 60 град широты до 70 град широты
c             qt=qt*(2.*nx-209.)/9.
c
c линейный рост в 3,5 раза от 60 град широты до 70 град широты
c             qt=qt*(5.*nx-509.)/36.
c
c             qt=qt*(9.*nx-476.)/64.
c             qt=qt*(9.*nx-663.)/48.
c             qt=qt*(9.*nx-803.)/34.
c             qt=qt*(9.*nx-476.)/128.
c           end if
c
c поток фотоэлектронов на 70 град широты увеличен в 6 раз
c           if(nx.eq.127)qt=qt*4.
c!!         if(nx.eq.127)qt=qt*6.
c!!         if(nx.eq.127)qt=qt*8.
c
c
c поток фотоэлектронов с 65 град широты увеличен в 4 раза
c           if(nx.ge.117)qt=qt*4.
c
c
c поток фотоэлектронов с 65 град широты увеличен в 4.5 раза
c           if(nx.ge.117)qt=qt*4.5
c
c
c поток фотоэлектронов с 65 град широты увеличен в 6 раз
c           if(nx.ge.117)qt=qt*6.0
c
c           qt=qt*(19.*nx+99.)/118.
c           qt=qt0*(1.-(1.+a)*exp(-a))*10.
c           a=(qo(1)+qo(nx))*.5
c           if(a.lt.10.)qt=qt*.1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            qt=qt/(snt*b1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           qt=qt/qomt
            qt=qt/qmax
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
   17 continue
      eht=1./dt
      do16m=i1,i2
        alt=ht(m)
        h=ch(m)
        vud=vdu(m)
        tet=tt(m)
        dv=dvn(alt,tet,vud)
        if(ne.eq.3)goto1
          tei=ti(m)
          ten=tn(m)
          voi=vio(m)
          if(alt.le.5.e8)o=co(m)
    1   continue
        if(ne.eq.1)goto2
          he=che(m)
    2   continue
        if(ne.eq.2)goto3
          if(alt.gt.1.e8)goto3
            n2=cn2(m)
            o2=co2(m)
    3   continue
        if(ne.eq.3.or.ne.eq.5)goto4
          vqn=vnq(m)
          vun=vnu(m)
          vvn=vnv(m)
          vvd=vdv(m)
    4   continue
        if(ne.eq.1.or.ne.eq.3)goto5
          oi=cio(m)
    5   continue
        if(ne.eq.2.or.ne.eq.3)goto6
          hi=cih(m)
    6   continue
        if(ne.eq.1.and.alt.le.1.e8.or.ne.eq.5.and.alt.le.1.e8)oq=qo(m)
        if(ne.eq.5.and.alt.le.1.e8)qmol=qsm(m)
        goto(7,7,8,13,13),ne
    7   continue
          if(alt.le.5.e8)la4=lamh(ten,tei)*o
c
          tv=ten
          if(nx.ge.45)then
            if(m.ge.1.and.m.le.23)tv=tvs(m)
            if(m.ge.nx-22.and.m.le.nx)tv=tvn(nx-m+1)
          end if
          la3=lamo(3,ten,vqn,vun,vvn,tei,voi,vud,vvd,tv)
c         la3=lamo(3,ten,vqn,vun,vvn,tei,voi,vud,vvd,te(m))
c         la3=lamo(3,ten,vqn,vun,vvn,tei,voi,vud,vvd)
c
          goto9
    8   continue
c         ga=n(he+)*eht
c         al=la5*n(o2)+la6*n(n2)+dv+eht,
          ga(m)=cihe1(m)*eht
c
c         ga(m)=0.
c
          if(dv.lt.0.)ga(m)=ga(m)-cihe1(m)*dv
          a=eht
          if(dv.ge.0.)a=dv+a
          if(alt.le.1.e8)a=a+la5*o2+la6*n2
          al(m)=a
          goto16
    9   continue
        goto(10,12),ne
   10   continue
c         ga=oq+la4*n(o)*n(h+)+n(o+)*eht
c         al=la1*n(o2)+la2*n(n2)+la3*n(h)+dv+eht,
          g=cio1(m)*eht
c
c         g=0.
c
          if(dv.lt.0.)g=g-cio1(m)*dv
          a=la3*h+eht
          if(dv.ge.0.)a=a+dv
          if(alt.gt.1.e8)goto11
c
          tv=ten
          if(nx.ge.45)then
            if(m.ge.1.and.m.le.23)tv=tvs(m)
            if(m.ge.nx-22.and.m.le.nx)tv=tvn(nx-m+1)
          end if
            la1=lamo(1,ten,vqn,vun,vvn,tei,voi,vud,vvd,tv)
            la2=lamo(2,ten,vqn,vun,vvn,tei,voi,vud,vvd,tv)
c           la1=lamo(1,ten,vqn,vun,vvn,tei,voi,vud,vvd,te(m))
c           la2=lamo(2,ten,vqn,vun,vvn,tei,voi,vud,vvd,te(m))
c           la1=lamo(1,ten,vqn,vun,vvn,tei,voi,vud,vvd)
c           la2=lamo(2,ten,vqn,vun,vvn,tei,voi,vud,vvd)
c
            g=g+oq+la5*o2*cihe(m)
            a=a+la1*o2+la2*n2
   11     continue
          if(alt.le.5.e8)g=g+la4*hi
          ga(m)=g
          al(m)=a
          goto16
   12   continue
c         ga=la3*n(h)*n(o+)+n(h+)*eht
c         al=la4*n(o)+dv+eht
          ga(m)=la3*h*oi+cih1(m)*eht
c
c         ga(m)=la3*h*oi
c
          if(dv.lt.0.)ga(m)=ga(m)-cih1(m)*dv
          a=eht
          if(dv.ge.0.)a=a+dv
          if(alt.le.5.e8)a=a+la4
          al(m)=a
          goto16
   13   continue
          mm=m-1
          mp=m+1
          if(m.eq.i1)mm=m
          if(m.eq.i2)mp=m
          hei=cihe(m)
          if(alt.le.1.e8)mi=cim(m)
          tee=te(m)
	!!!!!!!
	    if(tee.lt.ten) tee=ten	!!!!!
          pte=pite(oi,hi,hei,tee)
          ces=oi+hi+hei
          ce=ces
          if(alt.le.1.e8)ce=ce+mi
          bm=bdip(alt,tet)
          hp=ht(mp)
          hm=ht(mm)
          tp=tt(mp)
          tm=tt(mm)
          vhi=vih(m)
          vhei=vihe(m)
          dst=1./ds(hm,alt,hp,tm,tp)
          bmp=1./bdip(hp,tp)
          bmm=1./bdip(hm,tm)
          dvt=6.67e-1*dv
          w=oi*voi+hi*vhi+hei*vhei
          viop=vio(mp)
          viom=vio(mm)
          vihp=vih(mp)
          vihm=vih(mm)
          vihep=vihe(mp)
          vihem=vihe(mm)
          if(ne.eq.5)goto15
c     ga=(ptn*tn+ce*pte*te+pqj)/ces+ti*eht
c     al=(ptn+ce*pte+ra)/ces+dvt+eht
c     dvt=2/3*dv,
c     ces=n(o+)+n(h+)+n(he+),
c     ce=ces+n(m+)
            ptn=pitn(alt,o2,n2,o,h,he,ten,oi,hi,hei,tei)

            pqj=piqj(alt,o2,n2,o,h,he,ten,vqn,vun,vvn,
     *      oi,hi,hei,voi,vhi,vhei,tei,vud,vvd)
c           pqj=0.
            ces=1./ces
            b=ce*pte
            rao=oi*(viop*bmp-viom*bmm)*dst
            rah=hi*(vihp*bmp-vihm*bmm)*dst
            rahe=hei*(vihep*bmp-vihem*bmm)*dst
            ra=6.67e-1*bm*(rao+rah+rahe)
c           ra=0.
            g=(ptn*ten+pqj+b*tee)*ces+ti1(m)*eht
            if(ra.lt.0.)g=g-ti1(m)*ra*ces
            if(dvt.lt.0.)g=g-ti1(m)*dvt
c           a=(ptn+ra+b)*ces+eht
            a=(ptn+b)*ces+eht
            if(ra.ge.0.)a=a+ra*ces
            if(dvt.ge.0.)a=a+dvt
            ga(m)=g
            al(m)=a
            goto16
   15   continue
c     ga=(ptn*tn+pte*ti+ptd-ptd1-ptd2)+te*eht
c     al=ptn+pte+ra+dvt+eht
          q=0.
          if(alt.le.1.e8)q=oq
          ptn=pnte(alt,o2,n2,o,h,he,tee)
c старый вариант
c          if(alt.le.5.e8)ptn=ptn+pntrf(alt,o2,n2,o,ten,tee)
c старый вариант
c новый вариант
          if(alt.le.5.e8)ptn=ptn+pntrf1(alt,o2,n2,ten,tee)
c кроме Bailey, Balan
c          if(alt.le.5.e8)ptn=ptn+pntrf2(alt,o,ten,tee)
c кроме Bailey, Balan
c новый вариант
c
          tv=ten
          if(nx.ge.45)then
            if(m.ge.1.and.m.le.23)tv=tvs(m)
            if(m.ge.nx-22.and.m.le.nx)tv=tvn(nx-m+1)
          end if
c
          if(alt.le.1.e8)ptd1=petd12(o2,n2,ten,tee,tv)
c
c         if(alt.le.1.e8)ptd1=petd12(o2,n2,ten,tee)
c
          if(alt.le.5.e8)ptd2=petd3(o,ten,tee)
          ciom=cio(mm)
          cihm=cih(mm)
          cihem=cihe(mm)
          ciop=cio(mp)
          cihp=cih(mp)
          cihep=cihe(mp)
          cem=ciom+cihm+cihem
          cep=ciop+cihp+cihep
          if(hm.le.1.e8)cem=cem+cim(mm)
          if(hp.le.1.e8)cep=cep+cim(mp)
          vem=(ciom*viom+cihm*vihm+cihem*vihem)*bmm/cem
          vep=(ciop*viop+cihp*vihp+cihep*vihep)*bmp/cep
          ra=6.67e-1*bm*(vep-vem)*dst
c         ra=0.
          if(alt.le.1.e8)q=q+qmol
c          ptd=pgfel(q,ce,alt)
c
c дополнительный источник локального нагрева с 65 град широты
c         if(nx.ge.117)ptd=ptd*2.
c
c
c дополнительный источник локального нагрева с 70 град широты
ccc       if(nx.eq.127)ptd=ptd*2.3
c         if(nx.gt.127)ptd=ptd*3.
c
          ptd=pgfkr(q,ce,alt,o2,n2,o,h,he)
c          ptd=pgfkhaz(q,ce,alt,o2,n2,o,h,he)
c         if(alt.gt.7.e7)ptd=ptd+qt*bm*(qo(1)+qo(nx))*5.e-3/8.
          if(alt.gt.7.e7)ptd=ptd+qt*bm*(qo(iqo)+qo(nx-iqo+1))*.5
          g=ptn*ten+pte*tei+ptd+te1(m)*eht
c новый вариант
c Bailey, Balan
          g=g-pntrf2(alt,o,ten,tee)
c Bailey, Balan
c новый вариант
          if(ra.lt.0.)g=g-te1(m)*ra
          if(dvt.lt.0.)g=g-te1(m)*dvt
          if(alt.le.1.e8)g=g+ptd1
          if(alt.le.5.e8)g=g+ptd2
          ga(m)=g
c         a=ptn+pte+eht+ra
          a=ptn+pte+eht
          if(ra.ge.0.)a=a+ra
          if(dvt.ge.0.)a=a+dvt
          al(m)=a
   16 continue
     
      return
      end
!------------------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     N2 VIBRATIONAL QUANTA APPROXIMATION
c     Pavlov, A.V., The role of vibrationally excited nitrogen in the
c     formation of the mid-latitude negative ionospheric storms,
c     Annales Geophysicae, 1994, Vol.12, P. 554-564.
c
c     EQUATION (B.3)
c
c     INPUT:
c     h(nn)  - ALTITUDE(KM)
c
c     h(1)<140 km !!!  h(nn)>400 km !!!!  h(n)-h(n-1)<20 km !!!
c
c     yo(nn) - O NUMBER DENSITY(CM-3)
c     tn(nn) - TEMPERATURE
c     te(nn) - ELECTRON TEMPERATURE
c     cne(nn)- ELECTRON NUMBER DENSITY (CM-3)
c     nn - the number of altitude points
c
c    !!!!!  nn<200    !!!!!
c
c     OUTPUT:
c     alfa(nn) - VIBRATIONAL QUANTA
c     tv(nn)   - VIBRATIONAL TEMPERATURE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine altv(h,yo,tn,te,cne,alfa,tv,nn)
      dimension h(nn),yo(nn),tn(nn),te(nn),cne(nn),alfa(nn),
     *tv(nn),atv1(5),btv1(5),ctv1(5),atv2(5),btv2(5),ctv2(5)
      dimension clam(23),dl(23),w1(23),w2(23)
      double precision w(200),x1,w0,dn2,tdif,tvt,hn2,cl,cl1,w00
        data atv1/2.8,-2.745,-3.073,-4.,-4.469/
        data btv1/1.221e-3,4.454e-3,4.596e-3,4.994e-3,5.166e-3/
        data ctv1/-8.613e-8,-5.592e-7,-5.771e-7,-6.288e-7,-6.497e-7/
        data atv2/4.159,3.817,3.696,3.379,3.155/
        data btv2/7.042e-4,8.051e-4,8.303e-4,8.864e-4,9.232e-4/
        data ctv2/-4.166e-8,-4.814e-8,-4.959e-4,-5.301e-8,-5.505e-8/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     atv1,btv1,ctv1 for Te<4000 K,  atv2,btv2,ctv2 for Te>4000 K
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc  THE CALCULATION OF THE FREQUENCY AND LAMDA  cccc
        do 1 n=1,nn
        te0=te(n)
        w0=0.
        if(te0.lt.1200.) go to 4
        if(te0.gt.4000.) go to 3
        do 2 m=1,5
        x1=atv1(m)+btv1(m)*te0+ctv1(m)*te0*te0-16.d0
        x1=10.d0**x1
2       w0=w0+x1
        go to 4
3     continue
        do 5 m=1,5
        x1=atv2(m)+btv2(m)*te0+ctv2(m)*te0*te0-16.d0
        x1=10.d0**x1
5       w0=w0+x1
4     continue
c
       w(n)=w0*cne(n)
c
        g=981./(1.+h(n)/6378.)**2
        t=tn(n)
        hn2=826.e5*t/g/28.
        x1=yo(n)/1.d9
        dn2=9.69d7*t**0.724/x1
        tdif=hn2*hn2/dn2
        tvt=0.107*exp(-69.9/t**0.33)*x1
        tvt=1.d0/tvt
        clam(n)=2.d0*dsqrt(tdif/tvt)
        dl(n)=dsqrt(tdif*tvt)
1     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccc
         w0=0.
         w00=0.
         do 6 n=2,nn
         j=nn+1-n
         j1=j+1
        cl=clam(j)
        cl1=clam(j1)
        x1=(cl-cl1)/2.d0
        cl=dexp(-cl)
        cl1=dexp(-cl1)
        w0=w0+(w(j)*cl+w(j1)*cl1)*x1
        w00=w00+(w(j)/cl+w(j1)/cl1)*x1
         w1(j)=w0
         w2(j)=w00
6     continue
         w1(nn)=0.
         w2(nn)=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 7 n=1,nn
        x1=-3353./tn(n)
        teta=dexp(x1)
        cl=clam(n)
        cl1=dexp(cl)
        x1=1.d0/cl1
        alfa(n)=teta+dl(n)/cl*(w0*(cl1-x1)-cl1*w1(n)+x1*w2(n))
        cl=alfa(n)
	if(cl.le.0) then
	   print*,'altv n=',n,'cl=', cl
	   cl=10.d0
	end if 
        tv(n)=-3353./dlog(cl/(1.d0+cl))
7     continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!------------------------------------------------------------------------------
c cyclt2-vdrift-bdip ;
c        trubka-tube-bdip;
c                    lambet-bdip;
c                    vorvfv-bdip;
c                    backfv-bdip;
c                    ntnb-bdip;
c                    alga-bdip;
c
!      function bdip(alt,tet)
!      double precision d,t
!      data re/6371.02e5/,g10/-.30356/
!      t=tet
!      ro=1./(1.+alt/re)
!  900 format(' ',10g12.4)
!      roq=ro**3
!      d=dsqrt(1.d0+3.d0*dcos(t)**2)
!      bdip=-g10*roq*d
!      return
!      end
! Updated Mihail Melnik from Polar Geophysical Institute 14.02.2020
       function bdip(alt,tet)
        real,parameter :: re = 6371.02e5, g10 = .30356
        real ro,alt,tet
        real*8 t
        ro = re/(re+alt)
        t = dcos(dble(tet)) * dcos(dble(tet))
        bdip = g10 * ro * ro * ro * real(dsqrt(1.d0+3.d0*t ))
       end
!------------------------------------------------------------------------------
!      function ds(altm,alt,altp,tetm,tetp)
!      data re/6371.02e5/
!      ds=sqrt((altp-altm)**2+(alt+re)**2*(tetp-tetm)**2)
!      return
!      end
! Updated Mihail Melnik from Polar Geophysical Institute 14.02.2020
       function ds(altm,alt,altp,tetm,tetp)
       real,parameter :: re = 6371.02e5
        real altm,alt,altp,tetm,tetp
       ds=sqrt((altp-altm)*(altp-altm)+(alt+re)*(alt+re)*(tetp-tetm)*
     *  (tetp-tetm))
       end
!---------------------------------------------------------------------------
c       dvn :   lambet, alga
c
      function dvn(alt,tet,vdu)
      double precision s,c,cs,d,t
      data re/6371.02e5/
      t=tet
      s=dsin(t)
      c=dcos(t)
      cs=c*c
      d=dsqrt(1.d0+3.d0*cs)
      d=d**3
      dvn=-6.*s*(1.d0+cs)*vdu/((re+alt)*d)
      return
      end
!------------------------------------------------------------------------------
      subroutine ntnb(nx,ht,tt,cio,cih,cihe,cim,b1,snt,snb)

      dimension ht(*),tt(*),cio(*),cih(*),cihe(*),cim(*)
      is=18
      in=nx-18
      snb=0.
      snt=0.
      hm=ht(is)
      tm=tt(is)
      bm=bdip(hm,tm)
      b1=bm
      cem=cio(is)+cih(is)+cihe(is)
      if(hm.le.1.e8)cem=cem+cim(is)
      do 1 i=is,in
        ip=i+1
        hp=ht(ip)
        alt=(hm+hp)*.5
        tp=tt(ip)
        dst=ds(hm,alt,hp,tm,tp)*.5
        bp=bdip(hp,tp)
        cep=cio(ip)+cih(ip)+cihe(ip)
        if(hp.le.1.e8)cep=cep+cim(ip)
        snt=snt+(cem+cep)*dst
        snb=snb+(cem/bm+cep/bp)*dst
        hm=hp
        tm=tp
        bm=bp
        cem=cep
    1 continue
      snb=snb*b1
      return
      end
!------------------------------------------------------------------------------
      real function lamh(tn,ti)
      data c1/3.8e-11/
      a=sqrt(ti+tn*6.25e-2)
      lamh=c1*a
      return
      end
!------------------------------------------------------------------------------
      real function lamo(ne,tn,vnq,vnu,vnv,ti,vio,
     *vdu,vdv,te)
c    *vdu,vdv)
      data c11/2.82e-11/,c12/7.74e-12/,c13/1.073e-12/,c14/5.17e-14/,
     *c15/9.65e-16/,c21/1.533e-12/,c22/5.92e-13/,c23/8.6e-14/,
     *c24/2.73e-12/,c25/1.155e-12/,c26/1.483e-13/,c31/4.3e-11/,
     *c1/6.67e-1/,c2/4.28e-8/,c3/6.36e-1/,c4/4.08e-8/
      goto(1,1,4),ne
    1 continue
        vs=(vio-vnq)**2+(vdu-vnu)**2+(vdv-vnv)**2
        goto(2,3),ne
    2   continue
          t=c2*vs+c1*(ti-tn)+tn
          a=t/300.
          as=a*a
          aq=as*a
          af=aq*a
          lamo=c11-c12*a+c13*as-c14*aq+c15*af
          return
    3   continue
          t=c4*vs+c3*(ti-tn)+tn
c
          tv=te
c
c         tv=tn
c         tv=tn*1.15
cb        tv=tn*1.4
c         tv=tn*1.5
c         tv=tn*1.25
c         tv=tn*1.75
          a=t/300.
          as=a*a
c         if(t.le.1700.)lamo=c21-c22*a+c23*as
c         if(t.gt.1700.)lamo=c24-c25*a+c26*as
          if(t.le.1700.)bet0=c21-c22*a+c23*as
          if(t.gt.1700.)bet0=c24-c25*a+c26*as
          call betnv(bet,bet0,tv,t)
          lamo=bet
          return
    4 continue
        a=sqrt(tn+ti*6.25e-2)
        lamo=c31*a
        return
      end
!------------------------------------------------------------------------------
      subroutine betnv(bet,bet0,tv,tef)
c . . . calculation vibrational T N2
      dimension a(8),b(8)
      data a/3.39e-15,2.33e-14,3.02e-14,-2.74e-14,
     *      -3.84e-15,1.6e-14,-2.3e-14,  2.77e-14/
     *    ,b/3.72e-13,3.09e-11,1.92e-10,2.9e-10,
     *       5.85e-11,1.59e-10,1.19e-10,1.36e-10/
     *    ,e1/3353./
      tev=e1/tv
      x0=1.-exp(-tev)
      ski=0.
      do 1 i=1,8
       x=x0*exp(-tev*i)
c    . . . quasitrinor distribution
c      x=x0*exp(-i*tev+i*(i-1)*20.6/tef)
       ak=a(i)*tef+b(i)
       ski=ski+x*ak
    1 continue
      bet=bet0*x0+ski
      return
      end
!------------------------------------------------------------------------------
      function pite(cio,cih,cihe,te)
      data c1/3.70e-3/,c2/5.92e-2/,c3/1.48e-2/,s/-1.5/
      pite=(c1*cio+c2*cih+c3*cihe)*te**s
      return
      end
!------------------------------------------------------------------------------
      function pitn(alt,co2,cn2,co,ch,che,tn,cio,
     *cih,cihe,ti)
      real nu
      data c11/4.44e-10/,c12/5.00e-10/,c13/5.6e-11/,c14/2.42e-10/,
     *c15/2.11e-10/,c21/1.95e-10/,c22/2.34e-10/,c23/2.79e-10/,
     *c25/4.26e-10/,c31/4.68e-10/,c32/5.47e-10/,c33/4.35e-10/,
     *c34/7.58e-10/,c24/3.35e-10/,c35/8.75e-11/
      a=0.
      b=0.
      c=0.
      if(alt.gt.1.e8)goto1
        a=c11*co2+c12*cn2
        b=c21*co2+c22*cn2
        c=c31*co2+c32*cn2
    1 continue
      d=0.
      e=0.
      f=0.
      if(alt.gt.5.e8)goto2
        d=c13*nu(1,ti,tn)*co
        e=c23*co
        f=c33*co
    2 continue
c
      pitn1=a+d+c14*ch+c15*che
c
      pitn2=b+e+c24*nu(2,ti,tn)*ch+c25*che
c
      pitn3=c+f+c34*ch+c35*nu(1,ti,tn)*che
      pitn=cio*pitn1+cih*pitn2+cihe*pitn3
      return
      end
!------------------------------------------------------------------------------
      function petd12(co2,cn2,tn,te,tv)
      data c11/5.76e-9/,c12/3.902e3/,c13/4.38e2/,c14/4.56e-4/,
     *c15/2400./,c16/700./,c17/-3000./,c21/2.31e-8/,c22/2000./,
     *c23/1.06e4/,c24/7.51e3/,c25/1.1e-3/,c26/1800./,c27/3300./,
     *c28/1.233/,c29/1000./,c30/2.056e-4/,c31/4000./
c     tv=te
cb    tv=tn*1.4
c     tv=tn*1.5
c     tv=tn*1.25
c     tv=tn
c     tv=tn*1.15
c     tv=tn*1.75
      a=(te-tn)/(te*tn)
      b=exp(c17*a)-1.
      c=c13*tanh(c14*(te-c15))
      d=exp((c12+c)*(te-c16)/(c16*te))
c
      ptd1=c11*co2*d*b
      b=te-c29
      c=te-c31
      d=te-c26
      f=c23+c24*tanh(c25*d)
      g=c27+c28*b-c30*b*c
      a=(te-tv)/(te*tv)
      b=exp(-g*a)-1.
      c=exp(f*(te-c22)/(c22*te))
c
      ptd2=c21*cn2*c*b
      petd12=ptd1+ptd2
      return
      end
!------------------------------------------------------------------------------
      function petd3(co,tn,te)
      data c1/1.21e-8/,c2/3000./,c3/-22713./,c4/2.4e4/,
     *c5/3.e-1/,c6/1500./,c7/1.947e-5/,c8/4000./
      a=(te-tn)/(te*tn)
      if(a.lt.0) a=0
      b=exp(c3*a)-1.
      c=(te-c2)/(c2*te)
      d=te-c6
      e=te-c8
      f=c4+c5*d-c7*d*e
      a=exp(f*c)
      petd3=c1*co*a*b
      return
      end
!------------------------------------------------------------------------------
      function pgfkr(q,ce,alt,o2,n2,o,h,he)
      real n2
      data c1/1.238e3/
      pgfkr=0.
      if(alt.gt.1.e8)goto1
        cn=o2+n2+o+h+he
        cn=ce/cn
        a=alog10(cn)
        if(a.gt.-3.)a=-3.
        if(a.lt.-8.)a=-8.
        pgfkr=c1/ce*(9.+a)**2*q
    1 continue
      return
      end
!------------------------------------------------------------------------------
      function piqj(alt,co2,cn2,co,ch,che,tn,vnq,vnu,vnv,
     *cio,cih,cihe,vio,vih,vihe,ti,vdu,vdv)
      real nu
      data c1/2.44e-8/,c2/3.17e-13/,c3/1.01e-8/
      f=nu(3,ti,tn)
      vs=(vdu-vnu)**2+(vdv-vnv)**2
c     piqj1=vs*cio*pqji(1,alt,co2,cn2,co,ch,che,tn,ti)
c     piqj2=vs*cih*pqji(2,alt,co2,cn2,co,ch,che,tn,ti)
c     piqj3=vs*cihe*pqji(3,alt,co2,cn2,co,ch,che,tn,ti)
c     piqj5=0.
c     piqj4=0.
c     piqj6=0.
      piqj1=((vio-vnq)**2+vs)*cio*pqji(1,alt,co2,cn2,co,ch,che,tn,ti)
      piqj2=((vih-vnq)**2+vs)*cih*pqji(2,alt,co2,cn2,co,ch,che,tn,ti)
      piqj3=((vihe-vnq)**2+vs)*cihe*pqji(3,alt,co2,cn2,co,ch,che,tn,ti)
      piqj5=cio*cihe*(vio-vihe)**2*c1*f
      piqj4=cio*cih*(vio-vih)**2*c2
      piqj6=cih*cihe*(vih-vihe)**2*c3*f
      piqj=piqj1+piqj2+piqj3+piqj4+piqj5+piqj6
      return
      end
!------------------------------------------------------------------------------
c   nu:  lambet-tube-trubka-cyclt2
c   nu:  rik-lambet-tube-trubka-cyclt2
c   nu:  pitn-alga-tube-trubka-cyclt2
c   nu:  pqji-piqj- alga-tube-trubka-cyclt2
c   nu:  piqj- alga-tube-trubka-cyclt2
c
      real function nu(ne,ti,tn)
      data s1/.37/,s2/.38/,s3/-1.5/
      goto(1,1,4),ne
    1 continue
        t=ti+tn
        goto(2,3),ne
    2   continue
c
          a=t**s1
          goto5
    3   continue
c
          a=t**s2
          goto5
    4   continue
c		
  	 
          a=ti**s3
    5   continue
        nu=a
        return
      end
!------------------------------------------------------------------------------
      function pnte(alt,co2,cn2,co,ch,che,te)
      data c1/1.92e-14/,c2/9.13e-16/,c3/6.24e-15/,c4/4.94e-12/,
     *c5/1.26e-13/,c6/1.21e-4/,c7/3.6e-2/,c8/1.35e-4/
      ts=sqrt(te)
      pnte2=0.
      pnte1=0.
      pnte3=0.
      if(alt.gt.1.e8)goto1
        pnte2=c2*cn2*(1.-c6*te)*te
        pnte1=c3*co2*(1.+c7*ts)*ts
    1 continue
      if(alt.le.5.e8)pnte3=c1*co*ts
      pnte4=c4*ch*(1.-c8*te)*ts
      pnte5=c5*che*ts
      pnte=pnte1+pnte2+pnte3+pnte4+pnte5
      return
      end
!------------------------------------------------------------------------------
      function pntrf1(alt,co2,cn2,tn,te)
      data c1/5.34e-10/,c2/2.71e-10/,c3/2.63e-8/,c4/7.e-5/
      pntrf2=0.
      pntrf3=0.
      if(alt.gt.1.e8)goto1
        ts=1./sqrt(te)
        pntrf2=c1*co2*ts
        pntrf3=c2*cn2*ts
    1 continue
      pntrf1=pntrf2+pntrf3
      return
      end
!------------------------------------------------------------------------------
      function pntrf2(alt,co,tn,te)
      dimension eps(3),a(3),b(3),c(3),e(3),dx(3),ex(3)
      data c1/5.34e-10/,c2/2.71e-10/,c3/2.63e-8/,c4/7.e-5/
      data eps/0.02,0.028,0.08/,a/7.883e-6,9.466e-6,1.037e-8/
      data b/1.021,0.8458,1.633/,c/1.009,0.9444,1.466/
      data e/228.,326.,98./
c Pavlov, Berrington
c      d=5+exp(-326.6/tn)+3.*exp(-227.7/tn)
c      f=8.132e-9/d
c      pntrf2=c3/3.4e-12*f*co/(te*tn)
c Pavlov, Berrington
c Dalgarno, Degges
c      pntrf2=c3*co*(1.-c4*te)/tn
c Dalgarno, Degges
c a la Bailey
c       pntrf2=0.5*c3*co*(1.-c4*te)/tn
c a la Bailey
c Bailey, Balan
       t1=tn
       t0=tn
       z=5.+3.*exp(-228/t1)+exp(-326/t0)
       dx(1)=exp(-228./t1)
       dx(2)=exp(-326./t0)
       dx(3)=exp(-326./t0)
       ex(1)=exp(-228./te)
       ex(2)=exp(-326./te)
       ex(3)=exp(-98./te-228./t1)
       s=0.
       do i=1,3
         f=(1.+b(i))*dx(i)+(e(i)/te+1.+b(i))*ex(i)
         f=5.91e-9*(te-tn)*f+eps(i)*(ex(i)-dx(i))
         s=s+a(i)*c(i)*te**(b(i)-0.5)*f
       end do
       pntrf2=c3/3.4e-12*8.629e-6*co/z*s
c Bailey, Balan
      return
      end
!------------------------------------------------------------------------------
      subroutine algat(ne,nx,m,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                    vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *                    ti,te,vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                    al,ga,qomt,qmax,iqo,mass,NV)

      dimension ht(NV),tt(NV),co2(NV),cn2(NV),co(NV),ch(NV),
     *          che(NV),tn(NV),vnq(NV),vnu(NV),vnv(NV),cim(NV),
     *          cio(NV),cih(NV),cihe(NV),vio(NV),vih(NV),vihe(NV),
     *          vdu(NV),vdv(NV),cio1(NV),cih1(NV),cihe1(NV),
     *          ti1(NV),te1(NV),qo(NV),qsm(NV),
     *          mass(30)
!	* ,al(NV),ga(NV)
      real la1,la2,la3,la4,n2,lamh,lamo,mi,la5,la6
      data la5/1.e-9/,la6/1.4e-9/,e0/10./

      eht=1./dt
      alt=ht(m)
      h=ch(m)
      vud=vdu(m)
      tet=tt(m)
      dv=dvn(alt,tet,vud)
      tei=ti
      ten=tn(m)
      voi=vio(m)
!	print*,'algat m=',m,nv,nx,i1,i2,co
!	pause

      o=co(m)
      he=che(m)
      n2=cn2(m)
      o2=co2(m)
      if(ne.eq.5)goto4
        vqn=vnq(m)
        vun=vnu(m)
        vvn=vnv(m)
        vvd=vdv(m)
    4 continue
      oi=cio(m)
      hi=cih(m)
      if(ne.eq.5)oq=qo(m)
      if(ne.eq.5)qmol=qsm(m)
      mm=m-1
      mp=m+1
      if(m.eq.i1)mm=m
      if(m.eq.i2)mp=m
      hei=cihe(m)
      mi=cim(m)
      tee=te
      if (tee.lt.ten) TEE=TEN !!!!!
      pte=pite(oi,hi,hei,tee)
      ces=oi+hi+hei
      ce=ces
      ce=ce+mi
      bm=bdip(alt,tet)
      hp=ht(mp)
      hm=ht(mm)
      tp=tt(mp)
      tm=tt(mm)
      vhi=vih(m)
      vhei=vihe(m)
      dst=1./ds(hm,alt,hp,tm,tp)
      bmp=1./bdip(hp,tp)
      bmm=1./bdip(hm,tm)
      dvt=6.67e-1*dv
      w=oi*voi+hi*vhi+hei*vhei
      viop=vio(mp)
      viom=vio(mm)
      vihp=vih(mp)
      vihm=vih(mm)
      vihep=vihe(mp)
      vihem=vihe(mm)
      if(ne.eq.5)goto15
        ptn=pitn(alt,o2,n2,o,h,he,ten,oi,hi,hei,tei)
        pqj=piqj(alt,o2,n2,o,h,he,ten,vqn,vun,vvn,
     *  oi,hi,hei,voi,vhi,vhei,tei,vud,vvd)
        ces=1./ces
        b=ce*pte
        rao=oi*(viop*bmp-viom*bmm)*dst
        rah=hi*(vihp*bmp-vihm*bmm)*dst
        rahe=hei*(vihep*bmp-vihem*bmm)*dst
        ra=6.67e-1*bm*(rao+rah+rahe)
        g=(ptn*ten+pqj+b*tee)*ces+ti1(m)*eht
        if(ra.lt.0.)g=g-ti1(m)*ra*ces
        if(dvt.lt.0.)g=g-ti1(m)*dvt
        a=(ptn+b)*ces+eht
        if(ra.ge.0.)a=a+ra*ces
        if(dvt.ge.0.)a=a+dvt
        ga=g
        al=a
        goto16
   15 continue
        q=oq
        ptn=pnte(alt,o2,n2,o,h,he,tee)
c старый вариант
c        ptn=ptn+pntrf(alt,o2,n2,o,ten,tee)
c старый вариант
c новый вариант
        ptn=ptn+pntrf1(alt,o2,n2,ten,tee)
c кроме Bailey, Balan
c        ptn=ptn+pntrf2(alt,o,ten,tee)
c кроме Bailey, Balan
c новый вариант
        ptd1=petd12(o2,n2,ten,tee,ten)
        ptd2=petd3(o,ten,tee)
        ciom=cio(mm)
        cihm=cih(mm)
        cihem=cihe(mm)
        ciop=cio(mp)
        cihp=cih(mp)
        cihep=cihe(mp)
        cem=ciom+cihm+cihem
        cep=ciop+cihp+cihep
        cem=cem+cim(mm)
        cep=cep+cim(mp)
        vem=(ciom*viom+cihm*vihm+cihem*vihem)*bmm/cem
        vep=(ciop*viop+cihp*vihp+cihep*vihep)*bmp/cep
        ra=6.67e-1*bm*(vep-vem)*dst
        q=q+qmol
c        ptd=pgfel(q,ce,alt)
c
c дополнительный источник локального нагрева с 65 град широты
c       if(nx.ge.117)ptd=ptd*2.
c
c
c дополнительный источник локального нагрева с 70 град широты
ccc     if(nx.eq.127)ptd=ptd*2.3
c       if(nx.gt.127)ptd=ptd*3.
c
c        ptd=pgfkhaz(q,ce,alt,o2,n2,o,h,he)
        ptd=pgfkr(q,ce,alt,o2,n2,o,h,he)
        g=ptn*ten+pte*tei+ptd+te1(m)*eht
c новый вариант
c Bailey, Balan
        g=g-pntrf2(alt,o,ten,tee)
c Bailey, Balan
c новый вариант
        if(ra.lt.0.)g=g-te1(m)*ra
        if(dvt.lt.0.)g=g-te1(m)*dvt
        g=g+ptd1
        g=g+ptd2
        ga=g
        a=ptn+pte+eht
        if(ra.ge.0.)a=a+ra
        if(dvt.ge.0.)a=a+dvt
        al=a
   16 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine backfv(ne,i1,i2,ht,tt,cim,cio,cih,cihe,ti,te,
     *                  bet,gam,ga,al,lam)
      real lam
      dimension ht(*),tt(*),cim(*),cio(*),cih(*),ti(*),
     *          cihe(*),bet(*),gam(*),ga(*),al(*),lam(*),te(*)
      data pi/3.1415926/
      i3=i2-1
c     print*,' subrotine backfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
      do4mn=i1,i3
        m=i3-mn+i1
        mp=m+1
        mm=m-1
        if(m.eq.i1)mm=m
        alt=ht(m)
        bm=bdip(alt,tt(m))
        if(ne.gt.3)ce=cio(m)+cih(m)+cihe(m)
        if(ne.eq.5.and.alt.le.1.e8)ce=ce+cim(m)
        if(ne.gt.3)bm=bm/ce
        dst=ds(ht(mm),alt,ht(mp),tt(mm),tt(mp))
        dst=dst/(bm+bm)
        fl=lam(mp)-dst*ga(m)
c       if(ne.gt.3.and.mn.eq.i3)fl=0.
        fu=gam(m)-bet(m)*fl
c       if(m.ne.i1)lam(m)=fl+dst*al(m)*fu
        if(m.ne.i1)lam(m)=fl*(1.-bet(m)*dst*al(m))+dst*al(m)*gam(m)
        goto(1,2,3,5,6),ne
    1   continue
          cio(m)=fu
c         if(cio(m).lt.1.e-3)cio(m)=1.e-3
          goto4
    2   continue
          cih(m)=fu
c         if(cih(m).lt.1.e-3)cih(m)=1.e-3
          goto4
    3   continue
          cihe(m)=fu
c         if(cihe(m).lt.1.e-3)cihe(m)=1.e-3
          goto4
    5   continue
          ti(m)=fu
          goto4
    6   continue
          te(m)=fu
    4 continue
c     print*,' subrotine backfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
c     do i=i1,i2
c       write(*,10)ht(i)*1.e-5,tt(i)*180./pi,cio(i),cih(i),lam(i)
c  10   format(f8.0,f8.2,3(1pe12.2))
c     end do
c     print*,'backfv output'
      return
      end
!------------------------------------------------------------------------------
      subroutine backt(i1,i2,delu,bet,gam)

      dimension delu(*),bet(*),gam(*)
      i4=i2-1
      do1i=i1,i4
        m=i4-i+i1
        mp=m+1
        delu(m)=bet(m)*delu(mp)+gam(m)
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine forvfv(ne,i1,i2,ht,tt,cim,cio,cih,cihe,
     *alf,beta,bet,gam,ga,al,lam,vio,vih,vihe,col)
      real lam,lamp
      dimension ht(*),tt(*),cim(*),cio(*),cih(*),
     *cihe(*),bet(*),beta(*),gam(*),ga(*),al(*),
     *lam(*),vio(*),vih(*),vihe(*),col(*)
      data pi/3.14159265/,sigm/.93/
c     data pi/3.14159265/,sigm/.87/
c     print*,' subrotine forvfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
      i3=i1+1
      altm=ht(i1)
      altp=ht(i3)
      tetm=tt(i1)
      tetp=tt(i3)
      dstm=ds(altm,altm,altp,tetm,tetp)
      bm=bdip(altm,tetm)
      if(ne.gt.3)ce=cio(i1)+cih(i1)+cihe(i1)
      if(ne.eq.5.and.altm.le.1.e8)ce=ce+cim(i1)
      if(ne.gt.3)bm=bm/ce
      dstm=dstm/(bm+bm)
      i4=i2-1
      do1m=i1,i4
        g=ga(m)
        mp=m+1
        mm=m-1
        mn=mp+1
        if(m.eq.i1)mm=m
        if(m.eq.i4)mn=mp
        lamp=(lam(m)+lam(mp))*.5
cc      coll=(col(m)+col(mp))*.5
        alt=(altm+altp)*.5
        dst=1./ds(altm,alt,altp,tetm,tetp)
        alt=lamp*dst
        betp=(beta(m)+beta(mp))*.5
        ab=abs(betp)
        bp=(betp+ab)*.5
        lamp=alt+bp
        bm=(betp-ab)*.5
        betp=alt-bm
        dstp=ds(altm,altp,ht(mn),tetm,tt(mn))
        bgm=bdip(altm,tetm)
        bgmp=bdip(altp,tetp)
        bgmm=bdip(ht(mm),tt(mm))
        if(ne.gt.3)ce=cio(mp)+cih(mp)+cihe(mp)
        if(ne.eq.5.and.altp.le.1.e8)ce=ce+cim(mp)
        if(ne.gt.3)bgmp=bgmp/ce
        dstp=dstp/(bgmp+bgmp)
        a=al(mp)
        bmp=dstp
        bmm=dstm
        b=bet(m)
c       write(*,11)ht(m)*1.e-5,tt(m)*180./pi,lamp
c       write(*,11)ht(m)*1.e-5,tt(m)*180./pi,alt
c  11   format(f8.0,f8.2,1pe15.5)
        if(lamp.lt.1.)then
          alt=1.+b*betp
cc        alt=coll+b*betp
          dst=1./(lamp+bmp*alt*a)
          alf=dst*lamp
          bet(mp)=dst*alt
c         bmm=betp*(gam(m)+bmm*b*g/sigm)
c         goto(2,3,4),ne
c   2     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp*b*(bmm*al(m)*cio(m)+
c    *        .5*(cio(mp)*vio(mp)/bgmp-cio(mm)*vio(mm)/bgmm))+
c    *        .5*(cio(m)*vio(m)/bgm+cio(mp)*vio(mp)/bgmp)+
c    *        lamp*cio(mp)-betp*cio(m))
c           goto5
c   3     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp*b*(bmm*al(m)*cih(m)+
c    *        .5*(cih(mp)*vih(mp)/bgmp-cih(mm)*vih(mm)/bgmm))+
c    *        .5*(cih(m)*vih(m)/bgm+cih(mp)*vih(mp)/bgmp)+
c    *        lamp*cih(mp)-betp*cih(m))
c           goto5
c   4     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp*b*(bmm*al(m)*cihe(m)+
c    *        .5*(cihe(mp)*vihe(mp)/bgmp-cihe(mm)*vihe(mm)/bgmm))+
c    *        .5*(cihe(m)*vihe(m)/bgm+cihe(mp)*vihe(mp)/bgmp)+
c    *        lamp*cihe(mp)-betp*cihe(m))
c   5     continue
cc        gam(mp)=dst*bmm
          gam(mp)=dst*betp*(gam(m)+bmm*b*g)
        else
          alt=1./lamp+betp/lamp*b
cc        alt=coll/lamp+betp/lamp*b
          dst=1./(1.+bmp*alt*a)
          alf=dst
          bet(mp)=dst*alt
cc        bmm=betp/lamp*(gam(m)+bmm*b*g/sigm)
c         goto(6,7,8),ne
c   6     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp/lamp*b*(bmm*al(m)*cio(m)+
c    *        .5*(cio(mp)*vio(mp)/bgmp-cio(mm)*vio(mm)/bgmm))+
c    *        .5*(cio(m)*vio(m)/bgm+cio(mp)*vio(mp)/bgmp)/lamp+
c    *        cio(mp)-betp/lamp*cio(m))
c           goto9
c   7     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp/lamp*b*(bmm*al(m)*cih(m)+
c    *        .5*(cih(mp)*vih(mp)/bgmp-cih(mm)*vih(mm)/bgmm))+
c    *        .5*(cih(m)*vih(m)/bgm+cih(mp)*vih(mp)/bgmp)/lamp+
c    *        cih(mp)-betp/lamp*cih(m))
c           goto9
c   8     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp/lamp*b*(bmm*al(m)*cihe(m)+
c    *        .5*(cihe(mp)*vihe(mp)/bgmp-cihe(mm)*vihe(mm)/bgmm))+
c    *        .5*(cihe(m)*vihe(m)/bgm+cihe(mp)*vihe(mp)/bgmp)/lamp+
c    *        cihe(mp)-betp/lamp*cihe(m))
c   9     continue
cc        gam(mp)=dst*bmm
          gam(mp)=dst*betp/lamp*(gam(m)+bmm*b*g)
        end if
        altm=altp
        altp=ht(mn)
        tetm=tetp
        tetp=tt(mn)
        dstm=dstp
    1 continue
c     print*,' subrotine forvfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
c     print*,' alf=',alf
c     do i=i1,i2
c       write(*,10)ht(i)*1.e-5,tt(i)*180./pi,bet(i),gam(i)
c  10   format(f8.0,f8.2,2(1pe15.5))
c     end do
      return
      end
!------------------------------------------------------------------------------
      subroutine forvt(ne,i1,i2,ht,tt,cim,cio,cih,cihe,beta,bet,
     *                 gamm,ga,al,lam,betap,gap,
     *                 alp,lamp,betam,gam,alm,lamm,temp)
      real lam,lamp,lamm

      dimension ht(*),tt(*),cim(*),cio(*),cih(*),cihe(*),
     *          beta(*),bet(*),gamm(*),ga(*),al(*),lam(*),
     *          betap(*),betam(*),gap(*),gam(*),alp(*),alm(*),
     *          lamp(*),lamm(*),temp(*)
c     data dte/.5/
      data dte/.1/
      i3=i1+1
      i4=i2-1
      do1m=i3,i4
c       dte=temp(m)*5.e-5
        mp=m+1
        mm=m-1
        dtem=dte*2.*temp(m)
        dtemp=dte*2.*temp(mp)
        dtemm=dte*2.*temp(mm)
        alt=ht(m)
        altp=ht(mp)
        altm=ht(mm)
        tet=tt(m)
        tetp=tt(mp)
        tetm=tt(mm)
        bg=bdip(alt,tet)
        ce=cio(m)+cih(m)+cihe(m)
        if(ne.eq.5.and.alt.le.1.e8)ce=ce+cim(m)
        dst=ds(altm,alt,altp,tetm,tetp)
        bg=bg/(ce*dst)
        hp=(altp+alt)*.5
        dst=ds(alt,hp,altp,tet,tetp)
        gp=dst
        bm=beta(m)
        bp=abs(bm)
        a=lam(m)
        d=-(bm-bp)*.5/dst
        f=d*(temp(mp)-temp(m))
        e=bg*(a+lam(mp))/dst
        b=e+d
        f=f+e*(temp(mp)-temp(m))
        hm=(altm+alt)*.5
        dst=ds(altm,hm,alt,tetm,tet)
        gm=dst
        d=(bm+bp)*.5/dst
        f=f-d*(temp(m)-temp(mm))
        e=bg*(a+lam(mm))/dst
        f=f-e*(temp(m)-temp(mm))
        a=e+d
        c=a+b+al(m)
        c=c+(alp(m)-alm(m))/dtem*temp(m)-(gap(m)-gam(m))/dtem
        c=c+bg*(lamp(m)-lamm(m))*(temp(m)-temp(mm))/(gm*dtem)
        c=c-bg*(lamp(m)-lamm(m))*(temp(mp)-temp(m))/(gp*dtem)
        a=a-bg*(lamp(mm)-lamm(mm))/dtemm*(temp(m)-temp(mm))/gm
        b=b+bg*(lamp(mp)-lamm(mp))/dtemp*(temp(mp)-temp(m))/gp
        hp=betap(m)
        hm=betam(m)
        d=abs(hp)
        e=abs(hm)
        c=c+(hp-d-hm+e)*.5/dtem*(temp(mp)-temp(m))/gp
        c=c+(hp+d-hm-e)*.5/dtem*(temp(m)-temp(mm))/gm
        f=f+ga(m)-al(m)*temp(m)
        c=1./(c-a*bet(mm))
        bet(m)=b*c
        gamm(m)=(f+a*gamm(mm))*c
    1 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine lambet(ne,i1,i2,det,ht,tt,co2,cn2,co,ch,che,
     *           tn,vnq,cim,cio,cih,cihe,vio,vih,vihe,ti,te,col,
     *           vdu,vdv,beta,lam,vio1,vih1,vihe1,al,ga,
     *           cio1,cih1,cihe1,dolm)

      real lam,n2,nu4,nu5,nu,la

      dimension ht(*),tt(*),co2(*),cn2(*),co(*),ch(*),che(*),
     *          tn(*),vnq(*),cim(*),cio(*),cih(*),col(*),
     *          cihe(*),vio(*),vih(*),vihe(*),ti(*),te(*),
     *          vdu(*),vdv(*),beta(*),lam(*),vio1(*),vih1(*),
     *          vihe1(*),al(*),ga(*),cio1(*),cih1(*),cihe1(*)
      data re/6371.02e5/
c     print*,'lambet inpout'
c     print*,'ne=',ne,'i1=',i1,'i2=',i2,'dolm=',dolm
      eht=1./det
      do 7 m=i1,i2
        alt=ht(m)
        tet=tt(m)
        tei=ti(m)
        ten=tn(m)
        if(alt.le.1.e8) o2=co2(m)
        if(alt.le.1.e8) n2=cn2(m)
        if(alt.le.5.e8) o=co(m)
        h=ch(m)
        he=che(m)
        oi=cio(m)
        hi=cih(m)
        hei=cihe(m)
        tee=te(m)
        bm=bdip(alt,tet)
        ce=oi+hi+hei
        ces=ce
        if(alt.le.1.e8)ce=ce+cim(m)
        goto(1,1,1,6,6),ne
    1   continue
c  Classic electron condactivity
          mm=m-1
          mp=m+1
          if(m.eq.i1)mm=m
          if(m.eq.i2)mp=m
          hp=ht(mp)
          hm=ht(mm)
          dst=ds(hm,alt,hp,tt(mm),tt(mp))
          dst=1./dst
c  Classic electron condactivity end
          vud=vdu(m)
          gs=gst(alt,tet,vdv(m),vud,dolm)
          ci=cos(tet)
          ci=sin(tet)/sqrt(1.+3.*ci*ci)
          r=re+alt
          gs=gs-(2.*vud*oml(dolm)+omr(tet,dolm)*ci*
     *    (2.*vdv(m)-r*omt(tet,dolm)))
          dti=(ti(mp)-ti(mm))*dst
c  Classic electron condactivity
          dte=(te(mp)-te(mm))*dst
c  Classic electron condactivity end
          dnm=0.
          if(alt.lt.1.e8.and.hp.le.1.e8.and.hm.le.1.e8)
     *    dnm=(cim(mp)-cim(mm))*dst
          if(alt.le.1.e8.and.hp.gt.1.e8)
     *    dnm=(cim(m)-cim(mm))*dst*2.
          if(alt.le.1.e8.and.hm.gt.1.e8)
     *    dnm=(cim(mp)-cim(m))*dst*2.
          dnv=dvn(alt,tet,vud)
cccc
          rin=0.
cccc      rin=eht
cccc
c         rin=dnv
c         rin=dnv+eht
cccc
cccc      if(dnv.ge.0.)rin=rin+dnv*.5
cccc
          goto(2,3,4),ne
    2     continue
            r=rik(1,tei,ten,o2,n2,o,h,he,alt)
!	print*,'lambet ',ne,i1,i2,tei,ten
            nu5=nu(3,tei,ten)
            s1=sik(1,oi,hi,hei,nu5,nu4)
            s2=sik(2,oi,hi,hei,nu5,nu4)
            dnk=(cih(mp)-cih(mm))*dst
            dnl=(cihe(mp)-cihe(mm))*dst
            fob=1./oi
c           if(vio(m).ge.0.)dv=(vio(m)-vio(mm))*dst*2.
c           if(vio(m).lt.0.)dv=(vio(mp)-vio(m))*dst*2.
c           rin=rin+dv
c           rin=rin+(ga(m)*fob-al(m))+(1.-cio1(m)*fob)*eht
cccc
cccc        gs=gs+vio1(m)*eht
cccc        if(dnv.lt.0.)gs=gs-.5*dnv*vio(m)
cccc
            goto5
    3     continue
c
            r=rik(2,tei,ten,o2,n2,o,h,he,alt)
            nu4=nu(3,tei,ten)
            s1=sik(3,oi,hi,hei,nu5,nu4)
            s2=sik(4,oi,hi,hei,nu5,nu4)
            dnk=(cio(mp)-cio(mm))*dst
            dnl=(cihe(mp)-cihe(mm))*dst
            fob=1./hi
c           if(vih(m).ge.0.)dv=(vih(m)-vih(mm))*dst*2.
c           if(vih(m).lt.0.)dv=(vih(mp)-vih(m))*dst*2.
c           rin=rin+dv
c           rin=rin+(ga(m)*fob-al(m))+(1.-cih1(m)*fob)*eht
cccc
cccc        gs=gs+vih1(m)*eht
cccc        if(dnv.lt.0.)gs=gs-.5*dnv*vih(m)
cccc
            goto5
    4     continue
c
            r=rik(3,tei,ten,o2,n2,o,h,he,alt)
            nu5=nu(3,tei,ten)
            nu4=nu5
            s1=sik(5,oi,hi,hei,nu5,nu4)
            s2=sik(6,oi,hi,hei,nu5,nu4)
            dnk=(cio(mp)-cio(mm))*dst
            dnl=(cih(mp)-cih(mm))*dst
            fob=1./hei
c           if(vihe(m).ge.0.)dv=(vihe(m)-vihe(mm))*dst*2.
c           if(vihe(m).lt.0.)dv=(vihe(mp)-vihe(m))*dst*2.
c           rin=rin+dv
c           rin=rin+(ga(m)*fob-al(m))+(1.-cihe1(m)*fob)*eht
cccc
cccc        gs=gs+vihe1(m)*eht
cccc        if(dnv.lt.0.)gs=gs-.5*dnv*vihe(m)
cccc
    5     continue
          rps=1./((r+s1+s2+rin)*bm)
c         rps=1./bm
c         col(m) =r+s1+s2+rin
          r=r*vnq(m)
          dn=dnk+dnl
          if(alt.le.1.e8)dn=dn+dnm
          dt=tee/ce*dn+dti+dte
          r=(r+gs)*rps
          dt=dt*rps
          beta(m)=be_lambet(ne,dt,r,rps,vio(m),vih(m),vihe(m),s1,s2)
    6   continue
        if(ne.gt.3)w=oi*vio(m)+hi*vih(m)+hei*vihe(m)
c       if(ne.gt.3)w=0.
        if(ne.eq.4)beta(m)=w/ces
        if(ne.eq.5)beta(m)=w/ce
c  Classic electron condactivity
        lam(m)=la(ne,o2,n2,o,h,he,oi,hi,hei,ce,bm,rps,tei,tee,alt)
c  Classic electron condactivity end

    7 continue
c     print*,'lambet output'
      return
      end
!------------------------------------------------------------------------------
      function gst(alt,tet,vdv,vdu,dolm)
      double precision c,cs,t
      data re/6371.02e5/,g0/980.665/,om/7.272205e-5/
      t=tet
      c=dcos(t)
      cs=c*c
      r=1./(re+alt)
      b=2.*vdv*omt(tet,dolm)
cccc  b=0.
      a=omr(tet,dolm)
      a=(a*a-om*om)/r
      d=1./(1.+3.*cs)
      ds=sqrt(d)
      d=vdv*vdv*.5+(1.+cs)*d*vdu*vdu
cccc  d=0.
c     gst=2.*c*ds*(g0*re*r*re-3.*d)*r
      gst=2.*c*ds*((g0*re*r*re-3.*d)*r+a+b)
      return
      end
!------------------------------------------------------------------------------
c     omr :     gst, lambet
c
      function omr(tet,dolm)
      data ct0/9.8027118e-1/,st0/1.9765732e-1/,om/7.272205e-5/,
     *pi/3.14159265359/
      fi=dolm/180.*pi
      omr=om*(cos(tet)*ct0-sin(tet)*st0*cos(fi))
      return
      end
!------------------------------------------------------------------------------
c   omt :  gst, lambet, 
c
      function omt(tet,dolm)
      data ct0/9.8027118e-1/,st0/1.9765732e-1/,om/7.272205e-5/,
     *pi/3.14159265359/
      fi=dolm/180.*pi
      omt=-om*(sin(tet)*ct0+cos(tet)*st0*cos(fi))
      return
      end
!------------------------------------------------------------------------------
      real function la(ne,co2,cn2,co,ch,che,cio,cih,cihe,
     *ce,bm,rps,ti,te,alt)
      real la1,la2,la3,la4,la5
      data c11/5.2e6/,c21/8.31e7/,c31/2.08e7/,c41/1.84e8/,
     *c42/9.18e7/,c43/3.72e8/,c51/5.94e9/,c52/7.08e-12/,c53/2.54e-13/,
     *c54/9.02e-13/,c55/1.09e-16/,c56/1.09e-11/,c57/1.76e-10/,
     *c58/2.4e-14/,c59/1.8e-11/,s/2.5/
      goto(1,2,3,4,5),ne
    1 continue
c
        t=ti+cio/ce*te
        la=c11*t*rps
        return
    2 continue
c
        t=ti+cih/ce*te
        la=c21*t*rps
        return
    3 continue
c
        t=ti+cihe/ce*te
        la=c31*t*rps
        return
    4 continue
c
        t=ti**s/bm
        c=(c41*cihe+c42*cio+c43*cih)/ce
        la=t*c
        return
    5 continue
c
        t=te**s*c51/bm
        la1=0.
        la2=0.
        la3=0.
        if(alt.gt.1.e8)goto6
          sq=sqrt(te)
          la1=(c52+c53*sq)*co2
          la2=(c54-c55*te)*sq*cn2
    6   continue
        if(alt.le.5.e8)la3=c56*co
        la4=(c57-c58*te)*ch
        la5=c59*che
        c=1.+te*te/ce*(la1+la2+la3+la4+la5)
c       c=1.+sqrt(2.)+te*te/ce*(la1+la2+la3+la4+la5)
        la=t/c
        return
      end
!------------------------------------------------------------------------------
      function oml(dolm)
      data st0/1.9765732e-1/,om/7.272205e-5/,
     *pi/3.14159265359/
      fi=dolm/180.*pi
      oml=om*st0*sin(fi)
      return
      end
!------------------------------------------------------------------------------
      function rik(ne,ti,tn,co2,cn2,co,ch,che,alt)
      real nu
      data c11/6.67e-10/,c12/6.87e-10/,c14/1.29e-10/,c15/1.32e-10/,
     *c21/3.21e-9/,c22/3.39e-9/,c23/2.37e-9/,c25/1.06e-9/,
     *c31/2.11e-9/,c32/2.19e-9/,c33/1.09e-9/,c34/4.74e-10/,
     *c24/3.35e-10/,c35/8.77e-11/,c13/5.59e-11/
      goto(1,2,3),ne
    1 continue
        a=0.
        if(alt.le.1.e8)a=c11*co2+c12*cn2
        b=0.
        if(alt.le.5.e8)b=nu(1,ti,tn)*co*c13
        rik=a+b+c14*ch+c15*che
        return
    2 continue
        a=0.
        if(alt.le.1.e8)a=c21*co2+c22*cn2
        b=0.
        if(alt.le.5.e8)b=c23*co
        rik=a+b+nu(2,ti,tn)*ch*c24+c25*che
        return
    3 continue
        a=0.
        if(alt.le.1.e8)a=c31*co2+c32*cn2
        b=0.
        if(alt.le.5.e8)b=c33*co
        rik=a+b+c34*ch+nu(1,ti,tn)*che*c35
        return
      end
!------------------------------------------------------------------------------
      function sik(ne,cio,cih,cihe,nu5,nu4)
      real nu4,nu5
      data c1/2.47e-6/,c2/3.95e-5/,c3/1.256/,c4/1.9e-1/,
     *c5/7.6e-1/,c6/3.14e-1/
      goto(1,2,3,4,5,6),ne
    1 continue
        sik=c1*cih
        return
    2 continue
        sik=c4*nu5*cihe
        return
    3 continue
        sik=c2*cio
        return
    4 continue
        sik=c3*nu4*cihe
        return
    5 continue
        sik=c5*nu5*cio
        return
    6 continue
        sik=c6*nu4*cih
        return
      end
!------------------------------------------------------------------------------
      function be_lambet(ne,dt,r,rps,vio,vih,vihe,s1,s2)
      dimension c(3)
      data c/5.2e6,8.31e7,2.08e7/
      goto(1,2,3),ne
    1 continue
c
        a=(s1*vih+s2*vihe)*rps
        goto4
    2 continue
c
        a=(s1*vio+s2*vihe)*rps
        goto4
    3 continue
c
        a=(s1*vio+s2*vih)*rps
    4 continue
      be_lambet=c(ne)*dt-r-a
      return
      end
!------------------------------------------------------------------------------
      function pqji(ne,alt,co2,cn2,co,ch,che,tn,ti)
      real nu
      data c11/5.7e-17/,c12/5.61e-17/,c13/3.58e-18/,c14/9.72e-19/,
     *c15/3.39e-18/,c21/2.50e-17/,c22/2.62e-17/,c23/1.79e-17/,
     *c24/1.34e-18/,c25/6.83e-18/,c31/6.01e-17/,c32/6.14e-17/,
     *c33/2.79e-17/,c34/3.04e-18/,c35/1.41e-18/
      goto(1,2,3),ne
    1 continue
        a=0.
        b=0.
        if(alt.le.1.e8)a=c11*co2+c12*cn2
        if(alt.le.5.e8)b=c13*nu(1,ti,tn)*co
        pqji=a+b+c14*ch+c15*che
        return
    2 continue
        a=0.
        b=0.
        if(alt.le.1.e8)a=c21*co2+c22*cn2
        if(alt.le.5.e8)b=c23*co
        pqji=a+b+c24*nu(2,ti,tn)*ch+c25*che
        return
    3 continue
        a=0.
        b=0.
        if(alt.le.1.e8)a=c31*co2+c32*cn2
        if(alt.le.5.e8)b=c33*co
        pqji=a+b+c34*ch+c35*nu(1,ti,tn)*che
        return
      end