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
