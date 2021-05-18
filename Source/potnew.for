! subroutine potnew, conduc, degaps, facef, integr, intpar, jet_ph, jeth,
!      libman, psieq, wwpgl
      subroutine potnew(del0,ut,pkp,bmpy,bmpz,vsol,csol,na,
     *        ap,dst,al,ae,au,alt,nh,ddolgs,ddolgt,dtets,dtet,
     *        isp,readfl,kpar,kdf,kdu,nr,ldor,nfile,md,nl,ncs,ncn,
     *        par,pole,pot0,mast,tet1,tet2,tet3,pdpc,eps0,om,
     *        fac1,fac2,fac3,pglo,ids,nzapjet,lzapjet)
!  na=ntr
!  nl=idt
!  ncn=nl2
!  kparn=kpart0 ?
      dimension kdf(20),kdu(20),pole(ldor/4),alt(nh),mast(40),
     *          par(kpar,nh,ncs),pglo(kpar,nh,ncs,ids),pot0(na,nl,ncn),
     *          nzapjet(5),lzapjet(5)
      allocatable stt(:,:,:),sft(:,:,:),
     *   sff(:,:,:),alfa(:,:),psie(:),sitt(:,:),
     *   beta(:,:),gamma(:,:),delta(:,:),
     *   psi(:,:),eu(:,:),sid(:,:),sift(:,:),
     *   sfti(:,:),pst(:,:),psf(:,:),
     *   ev(:,:),siff(:,:),cur(:,:),
     *   pot(:,:,:),par1(:,:,:),curh(:,:,:)
      logical readfl
      data re/6371.02e5/,pi/3.14159265359/
      ! data cdu1/1./,cda1/1./,cdu2/1./,cda2/1./
      dpi=pi+pi
      readfl=.true.
      nfile=5
      md=1
      kparn=mast(16)
      allocate (stt(nl,na,ncn),sft(nl,na,ncn),
     *   sff(nl,na,ncn),alfa(nl,ncn),psie(nl),sitt(nl,ncn),
     *   beta(nl,ncn),gamma(nl,ncn),delta(nl,ncn),
     *   psi(nl,ncn),eu(nl,ncn),sid(nl,ncn),sift(nl,ncn),
     *   sfti(nl,ncn),pst(nl,ncn),psf(nl,ncn),
     *   ev(nl,ncn),siff(nl,ncn),cur(nl,ncn),
     *   pot(na,nl,ncn),par1(kparn,na,ncn),curh(nl,ncn,na))
      nvar=mast(17)
      ipr=mast(19)
      nvg=mast(27)
      nvpd=mast(28)
      nvfacf=mast(29)
      if(nvar.ne.1)then
        if(nvpd.eq.1)pdpc=10.06+14.44*pkp
      end if
      do 22 ll=1,nl
        dolg=(ll-1)*ddolgs
        call wwpgl (pglo,kpar,nh,ncs,ids,par,dolg,ddolgs)
c       call wws(readfl,nfile,kpar,dolg,ddolgs,tet,dtets,nh,kdf,
c    *  ldor,isp,md,par,pole,nr)
        call intpar(ncs,ncn,dtets,dtet,par,kpar,nh,par1,kparn,na)
        if(nvar.eq.1)goto100
          if(nvfacf.eq.0)goto100
            fi=dolg/180.*pi
            call magsm(ut,del0,fi,phism,1)
            tau=pi+phism
            if(tau.ge.dpi)tau=tau-dpi
            if(tau.lt.0.)tau=tau+dpi
            call facef(ncn,na,alt,dtet,par1,kparn,tet1,tet2,tet3,bmpy,
     *      tau,nvg,fac1,fac2,fac3,mast)
  100   continue
        call conduc(ll,ncn,nl,na,alt,nh,dtet,par1,stt,sft,sff,kparn)
        call integr(ll,ncn,nl,na,alt,nh,dtet,par1,stt,sft,sff,
     *  alfa,beta,sfti,pst,psf,kparn,gamma,delta,psi,
     *  sid,sift,siff,sitt)
        if(nvar.ne.2)then
          call psieq(ll,ncn,nl,na,alt,dtet,par1,stt,sft,sff,kparn,psi,
     *               psie)
        end if
   22 continue
      if(nvar.ne.2)then
        j1=(ncn+1)/2
        do ll=1,nl
	    lp=ll+1
	    lm=ll-1
	    if(ll.eq.1)lm=nl
	    if(ll.eq.nl)lp=1
	    psi(ll,j1)=psi(ll,j1)-(psie(lp)-psie(lm))/ddolgs*.5
	  end do
	end if 

      call degaps(nvar,ncn,nl,na,alt,nh,dtet,ddolgs,gamma,delta,psi,
     *pst,psf,sfti)
      do14k=1,na
        do14i=1,nl
          do14j=1,ncn
            pot(k,i,j)=0.
   14 continue
      if(nvar.eq.1)goto36
        a=(re+alt(16))**2
        j1=ncn/2
        do 135 j=1,j1
          j2=ncn-j+1
          tet=(j-1)*dtet
          if(tet.ne.tet3)goto51
            goto53
   51     continue
          if(tet.ne.tet2)goto52
            goto53
   52     continue
          if(tet.ne.tet1)go to 135
   53       continue
            e=tet/180.*pi
            ct=cos(e)
            st=sin(e)
            si=(ct+ct)/sqrt(1.+3.*ct*ct)
            b=-a*st*abs(si)
            do35i=1,nl
              fi=(i-1)*ddolgs/180.*pi
              call magsm(ut,del0,fi,phism,1)
              tau=pi+phism
              if(tau.ge.dpi)tau=tau-dpi
              if(tau.lt.0.)tau=tau+dpi

              if(tet.ne.tet1)goto34
                tau3=pi+pi
                tau1=mast(14)/12.*pi
                tau2=mast(15)/12.*pi
                if(tau.ge.0..and.tau.le.tau1)then
                  al=pi/(tau1-tau2)
                  be=pi*.5-al*tau1
                end if
                if(tau.gt.tau1.and.tau.le.tau2)then
                  al=pi/(tau2-tau1)
                  be=pi*.5-al*tau1
                end if
                if(tau.gt.tau2.and.tau.lt.tau3)then
                  al=pi/(tau1-tau2)
                  be=pi*1.5-al*tau2
                end if
                fc1=fac1*sin(al*tau+be)
                if(nvg.eq.1)goto76
                  psi(i,j)=psi(i,j)+fc1*b
                  psi(i,j2)=psi(i,j2)+fc1*b
                  goto35
   76           continue
                pot0(na,i,j)=pdpc*.5e11*sin(tau)
                pot0(na,i,j2)=pot0(na,i,j)
                pot(na,i,j)=pot0(na,i,j)
                pot(na,i,j2)=pot0(na,i,j)
                goto35
   34         continue
              if(tet.ne.tet2)goto54
                tau3=pi+pi
                tau1=mast(20)/12.*pi
                tau2=mast(21)/12.*pi
                if(tau.ge.0..and.tau.le.tau1)then
                  al=pi/(tau1-tau2)
                  be=pi*.5-al*tau1
                end if
                if(tau.gt.tau1.and.tau.le.tau2)then
                  al=pi/(tau2-tau1)
                  be=pi*.5-al*tau1
                end if
                if(tau.gt.tau2.and.tau.lt.tau3)then
                  al=pi/(tau1-tau2)
                  be=pi*1.5-al*tau2
                end if
                fc2=-fac2*sin(al*tau+be)
                psi(i,j)=psi(i,j)+fc2*b
                psi(i,j2)=psi(i,j2)+fc2*b
                goto35
   54         continue
              if(bmpy.eq.0.)goto35
                tau2=pi
                tau1=mast(22)/12.*pi
                tau3=mast(30)/12.*pi
                if(bmpy.lt.0.)goto55
                  if(tau.lt.tau1.or.tau.gt.tau3)goto35
                    if(tau.eq.tau2)goto35
                      if(tau.gt.tau2)goto56
   57                   continue
                        psi(i,j)=psi(i,j)-fac3*b
                        psi(i,j2)=psi(i,j2)+fac3*b
                        goto35
   56                 continue
                      psi(i,j)=psi(i,j)+fac3*b
                      psi(i,j2)=psi(i,j2)-fac3*b
                      goto35
   55           continue
                if(tau.lt.tau1.or.tau.gt.tau3)goto35
                  if(tau.eq.tau2)goto35
                    if(tau.gt.tau2)goto57
                      goto56
   35   continue
  135   continue
   36 continue
      it=0
   37 continue
      it=it+1
      call libman(ncn,nl,na,alt,nh,dtet,ddolgs,alfa,beta,gamma,
     *delta,psi,pot0,pot,om,tet1,tet2,tet3,nvg)
      eps=0.
      j1=ncn-1
      do40i=1,nl
        do40j=2,j1
          if(pot0(na,i,j).eq.0.)goto38
            e=(pot(na,i,j)-pot0(na,i,j))/pot0(na,i,j)
            e=abs(e)
            goto39
   38     continue
            e=1.
   39     continue
          if(e.gt.eps)eps=e
   40 continue
!!!!!!!!!!!!!!
      pot0=pot
!!!!!!!!!!!!!!
      if(eps.le.eps0)goto43
        goto 37
   43 continue
        readfl=.false.
        nanl=na*nl*ncn
        call wpotef(readfl,pot,nanl,kdf,kdu,ldor,isp)
	if(ipr.eq.1)then
      pmc=pi/180.
	j1=(ncn+1)/2
      call  jeth(pmc,j1,nl,ncn,dtet,pot0,ddolgs,alt,nh,na,nvar,
     *            sitt,siff,sift,sid,cur,nzapjet,lzapjet)
      call jet_ph(pmc,j1,nl,ddolgs,na,alt,ncn,dtet,nvar,
     *  pglo,kpar,nh,ncs,ids,sff,sft,pot0,curh,nzapjet,lzapjet)
	end if
      print 45,it
   45 format(' ',46x,i5,'-ая итерация')
      nc1=(ncn+1)/2
      pmaxn=pot(na,1,1)
      pminn=pot(na,1,1)
      pmaxs=pot(na,1,ncn)
      pmins=pot(na,1,ncn)
      do2j=1,nc1
        k=ncn-j+1
        do1i=1,nl
          pn=pot(na,i,j)
          ps=pot(na,i,k)
          if(pmaxn.lt.pn)pmaxn=pn
          if(pminn.gt.pn)pminn=pn
          if(pmaxs.lt.ps)pmaxs=ps
          if(pmins.gt.ps)pmins=ps
    1   continue
    2 continue
      dpn=(pmaxn-pminn)*1.e-11
      dps=(pmaxs-pmins)*1.e-11
      print  3,dps,dpn
    3 format(' ',35x,'dps=',1pe9.2,' kV',10x,'dpn=',1pe9.2,' kV')
c     call comef(nl,ncn,na,alt,nh,dtet,ddolgs,pot,eu,ev)
c     call peef(nl,ncn,na,dtet,ddolgs,pot,eu,ev)
      deallocate (stt,sft,
     *   sff,alfa,psie,sitt,
     *   beta,gamma,delta,
     *   psi,eu,sid,sift,
     *   sfti,pst,psf,
     *   ev,siff,cur,
     *   pot,par1,curh)

      return
      end
!------------------------------------------------------------------------------
      subroutine conduc(l,nc,nl,na,alt,nh,dtets,par,stt,sft,sff,kpar)
      dimension alt(nh),par(kpar,na,nc),stt(nl,na,nc),sft(nl,na,nc),
     *sff(nl,na,nc)
      data e/1.60219e-20/,oe/1.76e7/,oi/3.11e2/,g10/-.30356/,
     *re/6371.02e5/,ci1/4.23e-10/,ci2/4.28e-10/,ci3/2.58e-10/,
     *ce1/1.82e-10/,ce11/3.6e-2/,ce2/2.33e-11/,ce21/1.21e-4/,
     *ce3/2.8e-10/,pi/3.14159265359/
      ie=(nc+1)/2
      do2i=1,nc
        tet=(i-1)*dtets/180.*pi
        ct=cos(tet)
        st=sin(tet)
        sk=sqrt(1.+3.*ct*ct)
        si=(ct+ct)/sk
        ci=st/sk
        sis=si*si
        cis=ci*ci
        do2j=1,na
          roq=(re/(re+alt(j)))**3
          b=-g10*roq*sk
          ee=e/b
          ome=oe*b
          omi=oi*b
          cm=par(4,j,i)
          co2=par(1,j,i)
          cn2=par(2,j,i)
          co=par(3,j,i)
c         tn=par(5,j,i)
          te=par(5,j,i)
c         b=sqrt(tn)
          b=sqrt(te)
          fi=ci1*co2+ci2*cn2+ci3*co
c         fe=ce1*(1.+ce11*b)*b*co2+ce2*(1.-ce21*tn)*tn*cn2+ce3*b*co
          fe=ce1*(1.+ce11*b)*b*co2+ce2*(1.-ce21*te)*te*cn2+ce3*b*co
          bu=ome*ome
          cu=omi*omi
          co2=fi*fi
          cn2=fe*fe
          b=1./(bu+cn2)
          co=1./(cu+co2)
          sp=cm*ee*(omi*fi*co+ome*fe*b)
          sh=cm*ee*(cu*co-bu*b)
          s0=cm*ee*(omi/fi+ome/fe)
          if(i.ne.ie)goto1
            stt(l,j,i)=s0
            sft(l,j,i)=0.
            sff(l,j,i)=sp+sh*sh/sp
            goto2
    1     continue
          b=1./(s0*sis+sp*cis)
          stt(l,j,i)=s0*sp*b
          sft(l,j,i)=s0*sh*si*b
          sff(l,j,i)=(s0*sp*sis+(sp*sp+sh*sh)*cis)*b
    2 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine degaps(nvar,nc,nl,na,alt,nh,dtets,dfis,gamma,
     *delta,psi,pst,psf,sfti)
      dimension alt(nh),gamma(nl,nc),delta(nl,nc),psi(nl,nc),
     *pst(nl,nc),psf(nl,nc),sfti(nl,nc)
      data pi/3.14159265359/,re/6371.02e5/
      df=(dfis+dfis)/180.*pi
      dt=(dtets+dtets)/180.*pi
      df=1./df
      dt=1./dt
      r=re+alt(na)
      k=nc-1
      k1=(nc+1)/2
      do2i=2,k
        ip=i+1
        im=i-1
        do2j=1,nl
          jp=j+1
          jm=j-1
          if(j.eq.1)jm=nl
          if(j.eq.nl)jp=1
          gamma(j,i)=(sfti(jp,i)-sfti(jm,i))*df
          delta(j,i)=(sfti(j,ip)-sfti(j,im))*dt
          if(nvar.eq.2)goto1
            if(i.eq.k1)goto2
              psi(j,i)=r*((pst(j,ip)-pst(j,im))*dt+(psf(jp,i)-
     *        psf(jm,i))*df)
              goto2
    1     continue
          psi(j,i)=0.
    2 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine facef(nc,na,alt,dtets,par,kpar,tet1,tet2,tet3,bmpy,
     *           tau,nvg,fac1,fac2,fac3,mast)
      dimension alt(na),par(kpar,na,nc),mast(40)

	allocatable sp(:)
      data e /1.60219e-20/,oe/2.07e7/,oi/3.1e2/,g10/-.30356/,
     *     re/6371.02e5/,ci1/4.23e-10/,ci2/4.28e-10/,ci3/2.58e-10/,
     *     ce1/1.82e-10/,ce11/3.6e-2/,ce2/2.33e-11/,ce21/1.21e-4/,
     *     ce3/2.8e-10/,pi/3.14159265359/,
     *     alf0/4.2e-7/,t0/300./
	 
	allocate (sp(na))

      tet1s=180.-tet1
      tet2s=180.-tet2
      tet3s=180.-tet3
      alf=alf0*t0*2.
      do10i=1,nc
        tets=(i-1)*dtets
        if(tets.eq.tet1.or.tets.eq.tet1s)goto6
          if(tets.eq.tet2.or.tets.eq.tet2s)goto5
            if(tets.eq.tet3.or.tets.eq.tet3s)goto1
              goto10
    1       continue
            if(bmpy.eq.0.)goto10
c             tau1=10./12.*pi
              tau2=pi
c             tau3=14./12.*pi
              tau1=mast(22)/12.*pi
              tau3=mast(30)/12.*pi
              if(bmpy.lt.0.)goto4
                if(tau.lt.tau1.or.tau.gt.tau3)goto10
                  if(tau.eq.tau2)goto10
                    if(tau.gt.tau2)goto3
    2                 continue
                      fac=fac3
                      if(tets.eq.tet3)fac=-fac
                      goto7
    3               continue
                    fac=fac3
                    if(tets.eq.tet3s)fac=-fac
                    goto7
    4             continue
                  if(tau.lt.tau1.or.tau.gt.tau3)goto10
                    if(tau.eq.tau2)goto10
                      if(tau.gt.tau2)goto2
                        goto3
    5     continue
          fac=-fac2*sin(tau)
          goto7
    6   continue
        fac=fac1*sin(tau)
        if(nvg.eq.0)goto7
          fac=0.
    7   continue
        tet=tets/180.*pi
        ct=cos(tet)
        sk=sqrt(1.+3.*ct*ct)
        h=alt(1)
        sip=0.
        do8j=1,na
          roq=(re/(re+alt(j)))**3
          b=-g10*roq*sk
          ee=e/b
          ome=oe*b
          omi=oi*b
          cm=par(4,j,i)
          co2=par(1,j,i)
          cn2=par(2,j,i)
          co=par(3,j,i)
c         tn=par(5,j,i)
          te=par(5,j,i)
c         b=sqrt(tn)
          b=sqrt(te)
          fi=ci1*co2+ci2*cn2+ci3*co
c         fe=ce1*(1.+ce11* b)*b*co2+ce2*(1.-ce21*tn)*tn*cn2+ce3*b*co
          fe=ce1*(1.+ce11*b)*b*co2+ce2*(1.-ce21*te)*te*cn2+ce3*b*co
          bu=ome*ome
          cu=omi*omi
          co2=fi*fi
          cn2=fe*fe
          b=1./(bu+cn2)
          co=1./(cu+co2)
          sp(j)=cm*ee*(omi*fi*co+ome*fe*b)
          if(j.eq.1)goto8
            hp=alt(j)
            sip=(sp(j)+sp(j-1))*.5*(hp-h)
            h=hp
    8   continue
        do9j=1,na
          cm=par(4,j,i)
c         tn=par(5,j,i)
          te=par(5,j,i)
c         b=cm-sp(j)/(e*sip)*fac*tn/(alf*cm)
          b=cm-sp(j)/(e*sip)*fac*te/(alf*cm)
          if(b.lt.1.)b=1.
          par(4,j,i)=b
    9   continue
   10 continue
      deallocate (sp)
      return
      end
!------------------------------------------------------------------------------
      subroutine integr(l,nc,nl,na,alt,nh,dtets,par,stt,sft,sff,
     *           alfa,beta,sfti,pst,psf,kpar,gamma,delta,psi,
     *           sid,sift,siff,sitt)
      dimension alt(nh),par(kpar,na,nc),stt(nl,na,nc),sft(nl,na,nc),
     *          sff(nl,na,nc),alfa(nl,nc),beta(nl,nc),sfti(nl,nc),
     *          pst(nl,nc),psf(nl,nc),gamma(nl,nc),delta(nl,nc),
     *          psi(nl,nc),sid(nl,nc),siff(nl,nc),sift(nl,nc),
     *          sitt(nl,nc)
      data pi/3.14159265359/,re/6371.02e5/
      pmc=pi/180.
      do4i=1,nc
        tet=(i-1)*dtets*pmc
        ct=cos(tet)
        st=sin(tet)
        sk=sqrt(1.+3.*ct*ct)
        si=(ct+ct)/sk
        ci=st/sk
	  b=bdip(alt(1),tet)
        vt=par(7,1,i)
        vr=par(6,1,i)
        vf=par(8,1,i)
        tt=stt(l,1,i)
        ft=sft(l,1,i)
        h=alt(1)
        s1=0.
        s3=0.
        s5=0.
        s2=0.
        if(i.eq.1.or.i.eq.nc)goto1
          s4=0.
    1   continue
        ff=sff(l,1,i)
        do3j=2,na
          bu=bdip(alt(j),tet)
          vtu=par(7,j,i)
          vru=par(6,j,i)
          vfu=par(8,j,i)
          ttu=stt(l,j,i)
          ftu=sft(l,j,i)
          hu=alt(j)
          dh=(hu-h)*.5
          s1=s1+(tt+ttu)*dh
          c1=-(b*tt*vf+bu*ttu*vfu)*dh*si
          c2=-(b*ft*vt+bu*ftu*vtu)*dh*si
          c3=(b*ft*vr+bu*ftu*vru)*dh*ci
          s3=s3+c1+c2+c3
          s5=s5+(ft+ftu)*dh
          ffu=sff(l,j,i)
          s2=s2+(ff+ffu)*dh
          if(i.eq.1.or.i.eq.nc)goto2
            c4=-(b*ft*vf+bu*ftu*vfu)*dh*si
            c5=(b*ff*vt+bu*ffu*vtu)*dh*si
            c6=-(b*ff*vr+bu*ffu*vru)*dh*ci
            s4=s4+c4+c5+c6
    2     continue
          ff=ffu
          b=bu
          vt=vtu
          vr=vru
          vf=vfu
          tt=ttu
          ft=ftu
          h=hu
    3   continue
        gamma(l,i)=s1*1.e9
	  sitt(l,i)=s1
        alfa(l,i)=s1*st
        pst(l,i)=s3*st
        sfti(l,i)=s5
        delta(l,i)=s2*1.e9
	  sift(l,i)=s5
	  siff(l,i)=s2
        if(i.eq.1.or.i.eq.nc)goto4
          beta(l,i)=s2/st
          psf(l,i)=s4
	    sid(l,i)=s4
    4 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine intpar(ncs,nc,dtets,dtet,par,kpars,nh,par1,kpar,na)
      dimension par(kpars,nh,ncs),par1(kpar,na,nc)
      j1=1
      j2=2
      tet1=0.
      tet2=dtets
      j3=nc-1
      do4j=1,j3
        te=(j-1)*dtet
    1   continue
        if(te.ge.tet1.and.te.lt.tet2)goto2
          tet1=tet2
          tet2=tet2+dtets
          j1=j2
          j2=j2+1
          goto1
    2   continue
        dt=(te-tet1)/dtets
        do3ip=1,kpar
          if(ip.lt.4)ips=ip
c                PC
          if(ip.eq.4.or.ip.eq.5)ips=ip+2
          if(ip.gt.5)ips=ip+4
c                PC
c                Labtam
c         if(ip.eq.4)ips=ip+2
c         if(ip.ge.5)ips=ip+4
c                Labtam
          do3k=1,na
            df=par(ips,k,j2)-par(ips,k,j1)
            par1(ip,k,j)=par(ips,k,j1)+df*dt
    3   continue
    4 continue
      do5ip=1,kpar
        if(ip.lt.4)ips=ip
c                PC
        if(ip.eq.4.or.ip.eq.5)ips=ip+2
        if(ip.gt.5)ips=ip+4
c                PC
c                Labtam
c       if(ip.eq.4)ips=ip+2
c       if(ip.ge.5)ips=ip+4
c                Labtam
        do5k=1,na
          par1(ip,k,nc)=par(ips,k,ncs)
    5 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine jet_ph(pmc,j1,nl,ddolgs,na,alt,ncn,dtet,nvar,
     *  pglo,kpar,nh,ncs,ids,sff,sft,pot0,curh,nzapjet,lzapjet)

      dimension alt(nh),pglo(kpar,nh,ncs,ids),pot0(na,nl,ncn),
     *          curh(nl,ncn,na),nzapjet(5),lzapjet(5),sff(nl,na,ncn),
     *	  	  sft(nl,na,ncn) !(idt0,ntr0,nl20)
      data re/6371.02e5/
      open(16,file='jet_p_h',access='direct',recl=lzapjet(4))
      do i=1,nl
         dolg=(i-1)*ddolgs
         ip=i+1
         im=i-1
         if(i.eq.1)im=nl
         if(i.eq.nl)ip=1
         do k=1,na
            kp=k+1
            km=k-1
            if(k.eq.1)km=k
            if(k.eq.na)kp=k
	      h=alt(k)
            r=re+h
            do j=2,ncn-1
               tet=(j-1)*dtet
               t=tet*pmc
               ct=cos(t)
               st=sin(t)
               sk=sqrt(1.+3.*ct*ct)
               b=bdip(h,t)
               s=0.
               if(nvar.ne.2)then
                  si=2.*ct/sk
                  ci=st/sk
                  vnu=pglo(11,k,j,i)*si-pglo(10,k,j,i)*ci
                  vnv=pglo(12,k,j,i)*si
                  s=s+(sff(i,k,j)*vnu-sft(i,k,j)*vnv)*b
               end if
               efv=-(pot0(na,ip,j)-pot0(na,im,j))/(2.*ddolgs*pmc*r*st)
               if(j.ne.j1)then
                  efu=-(pot0(na,i,j+1)-pot0(na,i,j-1))/(2.*dtet*pmc*r)
               else
                  efu=-(pot0(na,i,j)-pot0(na,i,j-1))/(dtet*pmc*r)
               end if
               curh(i,j,k)=efv*sff(i,k,j)+efu*sft(i,k,j)+s
       nzapjet(4)=nzapjet(4)+1
                  write(16,rec=nzapjet(4))90.-tet,alt(k)*1.e-5,
     *  curh(i,j,k)*1.e11
	   end do
	end do
      end do
      close(16)
      return
      end
!------------------------------------------------------------------------------
      subroutine jeth(pmc,j1,nl,ncn,dtet,pot0,ddolgs,alt,nh,na,
     *             nvar,sitt,siff,sift,sid,cur,nzapjet,lzapjet)
      dimension nzapjet(5),lzapjet(5),cur(nl,ncn),pot0(na,nl,ncn),
     *           siff(nl,ncn),sift(nl,ncn),sid(nl,ncn),alt(nh),
     *           sitt(nl,ncn)
      data re/6371.02e5/
      open(17,file='jet_l_h',access='direct',recl=lzapjet(5))
	r=re+alt(16)
	do i=1,nl
	      ip=i+1
	      im=i-1
	      if(i.eq.1)im=nl
	      if(i.eq.nl)ip=1
	   do j=2,j1
	      js=ncn-j+1
	      tet=(j-1)*dtet*pmc
	      tets=(js-1)*dtet*pmc
	      st=sin(tet)
	      sts=sin(tets)
	      efv=-(pot0(na,ip,j)-pot0(na,im,j))/(2.*ddolgs*pmc*r*st)
	      if(j.ne.j1)then 
	        efu=-(pot0(na,i,j+1)-pot0(na,i,j-1))/(2.*dtet*pmc*r)
	        efus=-(pot0(na,i,js+1)-pot0(na,i,js-1))/(2.*dtet*pmc*r)
	        efvs=-(pot0(na,ip,js)-pot0(na,im,js))/
     *	                (2.*ddolgs*pmc*r*sts)
	      else 
		   efu=-(pot0(na,i,j)-pot0(na,i,j-1))/(r*dtet*pmc)
                    end if
		   cur(i,j)=efv*siff(i,j)+efu*sift(i,j)
		   if(j.ne.j1)cur(i,js)=efvs*siff(i,js)+efus*sift(i,js)
	       if(nvar.ne.2)then
	         cur(i,j)=cur(i,j)+sid(i,j)
	         if(j.ne.j1)cur(i,js)=cur(i,js)+sid(i,js)
	       end if
        end do
	end do
	jstart=2
 	jfinish=ncn-jstart+1	 
        do i=1,nl
	   phi=(i-1)*ddolgs
	   do j=jstart,jfinish
	     tet=(j-1)*dtet
             nzapjet(5)=nzapjet(5)+1
             write(17,rec=nzapjet(5))phi,90.-tet,cur(i,j)*1.e6,
     *       sitt(i,j)*1.e9,siff(i,j)*1.e9,sift(i,j)*1.e9
        end do
      end do
      close(17)
      return
      end
!------------------------------------------------------------------------------
      subroutine libman(nc,nl,na,alt,nh,dtets,dfis,alfa,
     *beta,gamma,delta,psi,pot0,pot,om,tet1,tet2,tet3,nvg)
      dimension alt(nh),alfa(nl,nc),beta(nl,nc),gamma(nl,nc),
     *delta(nl,nc),psi(nl,nc),pot0(na,nl,nc),pot(na,nl,nc)
      data pi/3.14159265359/,re/6371.02e5/,pmi0/14.9/
      df=dfis/180.*pi
      dt=dtets/180.*pi
      df=1./df
      dt=1./dt
      k=nc-1
      dfs=df*df
      dts=dt*dt
      df=df*.5
      dt=dt*.5
      m=(nc+1)/2
      do4i=1,nl
        ip=i+1
        im=i-1
        if(i.eq.1)im=nl
        if(i.eq.nl)ip=1
        do3j=2,k
          jp=j+1
          jm=j-1
          tet=(j-1)*dtets/180.*pi
          st=sin(tet)
          st=st*st
          pmi=(re+alt(na))/(re*st)
          p=0.
          f=0.
          if(pmi.ge.pmi0.or.pmi.lt.pmi0.and.j.eq.m)goto2
            l=nc-j+1
            if(j.gt.m)goto1
              lp=l+1
              lm=l-1
              e=alfa(i,l)
              g=(e+alfa(i,lp))*.5*dts
              h=gamma(i,l)*dt
              a=abs(h)
              b=g+h+a
              p=p+b
              f=f+b*pot0(na,i,lp)
              g=(e+alfa(i,lm))*.5*dts
              b=g-h+a
              p=p+b
              f=f+b*pot0(na,i,lm)
              e=beta(i,l)
              h=delta(i,l)*df
              g=(e+beta(ip,l))*.5*dfs
              a=abs(h)
              b=g-h+a
              p=p+b
              f=f+b*pot0(na,ip,l)
              g=(e+beta(im,l))*.5*dfs
              b=g+h+a
              p=p+b
              f=f+b*pot0(na,im,l)-psi(i,l)
              goto2
    1       continue
            pot(na,i,j)=pot(na,i,l)
            goto3
    2     continue
          if(j.ne.m)goto8
            ee=alfa(i,j)
            h=gamma(i,j)+(alfa(i,jp)-alfa(i,jm))*dt
            h=h/ee*dt
            a=abs(h)
            b=dts+h+a
            p=p+b
            f=f+b*pot0(na,i,jp)
            b=dts-h+a
            p=p+b
            f=f+b*pot0(na,i,jm)
            e=beta(i,j)
            h=delta(i,j)*df
            g=(e+beta(ip,j))*.5*dfs
            a=abs(h)
            b=(g-h+a)/ee
            p=p+b
            f=f+b*pot0(na,ip,j)
            g=(e+beta(im,j))*.5*dfs
            b=(g+h+a)/ee
            p=p+b
            f=f+b*pot0(na,im,j)
            f=f-psi(i,j)/ee
            goto9
    8     continue
          if(nvg.eq.0)goto7
            if(j.eq.4.or.j.eq.34)goto3
    7     continue
          e=alfa(i,j)
          h=gamma(i,j)*dt
          g=(e+alfa(i,jp))*.5*dts
          a=abs(h)
          b=g+h+a
          p=p+b
          f=f+b*pot0(na,i,jp)
          g=(e+alfa(i,jm))*.5*dts
          b=g-h+a
          p=p+b
          f=f+b*pot0(na,i,jm)
          e=beta(i,j)
          h=delta(i,j)*df
          g=(e+beta(ip,j))*.5*dfs
          a=abs(h)
          b=g-h+a
          p=p+b
          f=f+b*pot0(na,ip,j)
          g=(e+beta(im,j))*.5*dfs
          b=g+h+a
          p=p+b
          f=f+b*pot0(na,im,j)
          f=f-psi(i,j)
    9     continue
          pot(na,i,j)=om/p*f+(1.-om)*pot0(na,i,j)
    3   continue
    4 continue
      k=nc-1
      s1=0.
      s2=0.
      do5i=1,nl
        s1=s1+pot(na,i,2)
        s2=s2+pot(na,i,k)
    5 continue
      s1=s1/nl
      s2=s2/nl
      do6i=1,nl
        pot(na,i,1)=s1
        pot(na,i,nc)=s2
    6 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine psieq(l,nc,nl,na,alt,dtets,par,stt,sft,sff,kpar,psi,
     *              psie)
      dimension alt(na),par(kpar,na,nc),stt(nl,na,nc),psi(nl,nc),
     *         sft(nl,na,nc),sff(nl,na,nc),psie(nl)
      data g10/-.30356/,pi/3.14159265359/,re/6371.02e5/
      i=(nc+1)/2
c      roq=(re/(re+alt(1)))**3
c      b=-g10*roq
      b=bdip(alt(1),pi*.5)
      r=re+alt(16)
      vf=par(8,1,i)
      vr=par(6,1,i)
      tt=stt(l,1,i)
      ff=sff(l,1,i)
      h=alt(1)
      s1=0.
      s2=0.
c         substorm
c      nm=na-4
c         substorm
c         matias firster
      nm=na-5
c         matias firster
c	nm=na-4
      if(l.eq.1)then
        print2,' источник динамо-поля на экваторе',
     *         ' рассчитывается до ',nm,'-ой высотной точки'
    2   format(a34,a19,i3,a18)
      end if
      do1j=2,nm
c        roq=(re/(re+alt(j)))**3
c        bu=-g10*roq
        bu=bdip(alt(j),pi*.5)
        vfu=par(8,j,i)
	vru=par(6,j,i)
        ttu=stt(l,j,i)
	ffu=sff(l,j,i)
	hu=alt(j)
        dh=(hu-h)*.5
        s1=s1+(b*tt*vf+bu*ttu*vfu)*dh*2.
	s2=s2+(b*ff*vr+bu*ffu*vru)*dh
        b=bu
        vf=vfu
	vr=vru
        tt=ttu
	ff=ffu
        h=hu
    1 continue
      psi(l,i)=r*s1
c	psie(l)=0.
      psie(l)=r*s2
      return
      end
!------------------------------------------------------------------------------
c                WWPGL  (cyclt1:potnew:wwpgl )
c             -----------------------------------
      subroutine wwpgl (pglo,kpars,nh,its,ids,par,dolg,ddolgs)
c                vhod - pglo,kpars,nh,its,ids,par,dolg,ddolgs
c                vihod - par
      dimension pglo(kpars,nh,its,ids),par(kpars,nh,its)
      nn=dolg
      nd=ddolgs
      if(nn.ge.360)nn=nn-360
      i=nn/nd+1
      do 2 j=1,its
        do 3 k=1,nh
          do 4 n=1,kpars
             par(n,k,j)=pglo(n,k,j,i)
    4     continue
    3   continue
    2 continue
      return
      end
!------------------------------------------------------------------------------