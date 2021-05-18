! subroutine potnewmf2, cintn, conducn, jet_pB2, jetB2_80, jetB2_h, jetB2,
!      degapn, sftpoln, libmanmf, libme, libmsn, libmu
      subroutine potnewmf2(del0,ut,pkp,bmpy,bmpz,vsol,csol,na,
     *        ap,dst,al,ae,au,alt,nh,ddolgs,ddolgt,dtets,dtet,
     *        isp,readfl,kpar,kdf,kdu,nr,ldor,nfile,md,nl,ncs,ncn,
     *        par,pole,pot0,mast,tet1,tet2,tet3,pdpc,eps0,om,
     *        fac1,fac2,fac3,pglo,ids,nzapjet,lzapjet)
      dimension kdf(20),kdu(20),par(kpar,nh,ncs),pole(ldor/4),alt(nh),
     *          mast(40),pot0(na,nl,ncn),pglo(kpar,nh,ncs,ids),
     *          nzapjet(5),lzapjet(5)
      ! 19 - число силовых линий	(ncn+1)/2 - для нечетных
	! 16 - na число высотных точек до 175 км
	! 41 число точек по широте + 4 доп точки для экватора - только для 5 град
	! 672=na*(41+1)
	! 32 = 16*2 - 2 координаты для 16 точек (от 80 до 175)
	! 1216 = 16*2 * 19*2 точки координат для указанной сетки
	! 4864 = 1216*4 (на чило интерполяционных параметров для расчета проводимости)
	allocatable  pefgl(:),sip(:,:),sih(:,:),sib(:,:),
     *             alfa(:,:),beta(:,:),gamma(:,:),delta(:,:),
     *             psi(:,:),sfti(:,:),pst(:,:),psf(:,:),
     *             Eu(:,:),Ev(:,:),pef0(:,:),pef(:,:),
     *	         aeb(:),beb(:),geb(:),cur(:,:),
     *             park(:),u(:),pari(:),q(:,:),ntsl(:),
     *             rads(:),pot(:,:,:),sfti2(:,:),sfti3(:,:),
     *             curh(:,:,:),cur80(:,:)

	logical readfl
      data re/6371.02e5/,pi/3.14159265359/
	dpi=pi+pi 
      nc=ncn+4	      ! 37+4=41 
	n19=(ncn+1)/2     ! only  for odd point ncn
	nrpef= na* (nc+1) ! 672
      naD=na*2          ! 16*2=32
	nark=naD*n19*2	  ! 1216
      kparn=mast(16)	  ! 4 
	ni=kparn*nark	  ! 4864

	allocate (pefgl(nrpef),sip(nl,nc),sih(nl,nc),sib(nl,nc),
     *          alfa(nl,nc),beta(nl,nc),gamma(nl,nc),delta(nl,nc),
     *          psi(nl,nc),sfti(nl,nc),pst(nl,nc),psf(nl,nc),
     *          Eu(nl,ncn),Ev(nl,ncn),pef0(nl,nc),pef(nl,nc),
     *		  aeb(nl),beb(nl),geb(nl),cur(nl,ncn),
     *          park(nark),u(n19+1),pari(ni),q(naD,n19),ntsl(n19),
     *          rads(nh),pot(na,nl,ncn),sfti2(nl,nc),sfti3(nl,nc),
     *          curh(nl,ncn,na),cur80(nl,nc))

      !readfl=.true.
      !nfile=5
      !md=1
      
      nvar=mast(17)
      ipr=mast(19)
      nvg=mast(27)
      nvpd=mast(28)
      nvfacf=mast(29)
      if(nvar.eq.1)goto12
        if(nvpd.eq.1)pdpc=10.06+14.44*pkp
   12 continue
      call iseti(nh,na,nui,u,nxi,nqi,q,rads,ns,park,ntsl)
      do22ll=1,nl
        dolg=(ll-1)*ddolgs
        call wwpgl (pglo,kpar,nh,ncs,ids,par,dolg,ddolgs)
        call gu(ll,rads,nh,par,kpar,ncs,nl,nvar,aeb,beb,geb)
        call insti(ntsl,nqi,par,pari,ni,park,kparn,rads,nh,ns,ncs,
     *  dtets,kpar)
        call sftpoln(ll,rads,nh,na,par,kpar,ncs,nc,sfti,nl)
        call cintn(ll,nvar,nqi,ntsl,kparn,nc,q,nxi,park,ns,pari,ni,
     *  nl,alfa,beta,sfti,sfti2,sfti3,psf,pst,rads,nh,par,kpar,ncs,
     *  sip,sih,sib)
c       if(nvar.eq.1)goto100
c         if(nvfacf.eq.0)goto100
c           fi=dolg/180.*pi
c            call magsm(ut,del0,fi,phism,1)
c           tau=pi+phism
c           if(tau.ge.dpi)tau=tau-dpi
c           if(tau.lt.0.)tau=tau+dpi
c           call facef(ncn,na,alt,dtet,par1,kparn,tet1,tet2,tet3,bmpy,
c    *      tau,nvg,fac1,fac2,fac3)
c 100   continue
   22 continue
c      if(ipr.ne.1)goto49
c        call pecondn(nl,nc,dtet,ddolgs,sip,sih,alt,nh)
c   49 continue
       call degapn(nvar,nc,nl,dtet,ddolgs,sfti,sfti2,sfti3,
     *       psf,pst,gamma,delta,psi,u,nui,rads,nh,na)
      do14k=1,na
        do14i=1,nl
          do14j=1,ncn
            pot(k,i,j)=0.
   14 continue
       pef0=0.
       pef=0.
!      do32i=1,nl
!        do32j=1,nc
!          pef0(i,j)=0.
!          pef(i,j)=0.
!   32 continue
      if(nvar.eq.1)goto36
        a=(re+alt(16))**3/re
        j1=ncn/2
        do135j=1,j1
          j2=nc-j+1
          tet=(j-1)*dtet
          if(tet.ne.tet3)goto51
            goto53
   51     continue
          if(tet.ne.tet2)goto52
            goto53
   52     continue
          if(tet.ne.tet1)goto135
   53       continue
            e=tet/180.*pi
            ct=cos(e)
            st=sin(e)
            sk=sqrt(1.+3.*ct*ct)
            si=(ct+ct)/sk
            b=-a/sk
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
                pef0(i,j)=pdpc*.5e11*sin(tau)
                pef(i,j)=pef0(i,j)
                pef0(i,j2)=pef0(i,j)
                pef(i,j2)=pef0(i,j)
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
      
	! iteration cycle
	it=0

   37 continue
		it=it+1
		call libmanmf(nc,nl,na,alt,nh,dtet,ddolgs,alfa,beta,gamma,
     *       delta,psi,pef0,om,tet1,tet2,tet3,aeb,beb,geb,pef,u,nui,nvg)
		eps=0.
		j1=nc-1
		do  i=1,nl
		  do  j=2,j1
			if(pef0(i,j).eq.0.)goto38
			  e=(pef(i,j)-pef0(i,j))/pef0(i,j)
			  e=abs(e)
			  goto39
   38       continue
			  e=1.
   39       continue
			if(e.gt.eps)eps=e
	       end do
         end do
!		  do 41 i=1,nl
!			do 41 j=1,nc
!			  pef0(i,j)=pef(i,j)
!	   41   continue
                pef0=pef

		if(eps.le.eps0) goto 43
		  goto37
   43 continue
	ami=pef(1,1)
      ama=pef(1,1)
      do801i=1,nl
        do801j=2,j1
          if(pef(i,j).lt.ami)ami=pef(i,j)
          if(pef(i,j).gt.ama)ama=pef(i,j)
  801 continue
      ami=(ama+ami)*.5
      do802i=1,nl
        do802j=1,nc
          pef(i,j)=pef(i,j)-ami
  802 continue
      j1=(ncn+1)/2
      do600i=1,nl
        do600j=1,ncn
          if(j.le.j1)pot(na,i,j)=pef(i,j)
          if(j.gt.j1)pot(na,i,j)=pef(i,j+4)
          if(j.le.j1)pot0(na,i,j)=pef(i,j)
          if(j.gt.j1)pot0(na,i,j)=pef(i,j+4)
  600 continue
c     call sumef(nl,ncn,na,alt,nh,dtet,ddolgs,pot)
      do803i=1,nl
        do803j=1,nc
          pef0(i,j)=pef(i,j)
  803 continue
      j1=(ncn+1)/2
      j2=(nc+1)/2
      nma=na-1
      do i=1,nl
        pot0(1,i,j1)=pef(i,j2)
        do k=1,nma
          pot0(k,i,1)=pef(i,1)
          pot0(k,i,ncn)=pef(i,nc)
        end do
      end do
      do idol=1,nl
        do j=1,nqi
          call filpef(idol,j,ntsl,nqi,pef,nl,nc,pefgl,nrpef)
        end do
        call intspef(idol,pefgl,nrpef,rads,nh,park,ns,ntsl,nqi,
     *       dtet,u,nui,ddolgs,dtets,pot0,na,nl,ncn,nc)
      end do

c      do k=1,na
c         print*,' k=',k
c         print*,'   Long','  Lat','  Pot'
c         do j=1,nl
c           dolg=(j-1)*15.
c            do i=1,37
c               phi=90.-(i-1)*5.
c              print*,dolg,phi,pot0(k,j,i)*1.e-11
c            end do
c         end do
c     end do

      readfl=.false.
      nanl=na*nl*ncn
      call wpotef(readfl,pot0,nanl,kdf,kdu,ldor,isp)

      if(ipr.eq.1)then   ! ipr=mast(19)
         pmc=pi/180.
		r=re+alt(16)
		call jetB2(pmc,j1,r,nl,ncn,dtet,pot0,na,ddolgs,pef,
     *         nc,u,cur,sip,sih,sib,nzapjet,lzapjet)
	 		 
		call jetB2_80(pmc,j1,j2,r,alt,nh,nl,ncn,dtet,pot0,na,
     *         ddolgs,pef,nc,u,cur80,sip,sih,sib,nzapjet,lzapjet)
		call jetB2_h(pmc,j1,r,alt,nl,ncn,dtet,nvar,pot0,na,ddolgs,
     *         pglo,kpar,nh,ncs,ids,cur,nzapjet,lzapjet)
		call jet_pB2(pmc,j1,nl,ddolgs,na,alt,ncn,dtet,nvar,
     *         pglo,kpar,nh,ncs,ids,pot0,curh,nzapjet,lzapjet)
      end if
      print 45,it
   45 format(' potnewmf2 ',i7,' iteration')
      nc1=(ncn+1)/2
      pmaxn=pot0(na,1,1)
      pminn=pot0(na,1,1)
      pmaxs=pot0(na,1,ncn)
      pmins=pot0(na,1,ncn)
      do2j=1,nc1
        k=ncn-j+1
        do1i=1,nl
          pn=pot0(na,i,j)
          ps=pot0(na,i,k)
      !    if(ISNAN(ps).or.ISNAN(ps)) then
      !      print*,'ps=',ps,' pn=',pn,i,j,k
      !    end if
          if(pmaxn.lt.pn)pmaxn=pn
          if(pminn.gt.pn)pminn=pn
          if(pmaxs.lt.ps)pmaxs=ps
          if(pmins.gt.ps)pmins=ps
    1   continue
    2 continue
      dpn=(pmaxn-pminn)*1.e-11
      dps=(pmaxs-pmins)*1.e-11
      print  3,dps,dpn
 3    format(' ','dps=',1pe9.2,' kV',' dpn=',1pe9.2,' kV')
      
      deallocate (pefgl,sip,sih,sib,
     *           alfa,beta,gamma,delta,
     *           psi,sfti,pst,psf,
     *           Eu,Ev,pef0,pef,
     *	       aeb,beb,geb,cur,
     *           park,u ,pari ,q ,ntsl ,
     *           rads,pot,sfti2,sfti3,
     *           curh,cur80)

      return
      end
!------------------------------------------------------------------------------
      subroutine cintn(kk,nvar,nqi,ntsl,kparn,ncn,qq,nxi,park,ns,
     *                 pari,ni,nl,alfa,beta,sfti1,sfti2,sfti3,
     *                 psf,pst,rads,nh,par,kpars,ncs,sip,sih,sib)
      dimension msum(45),ntsl(nqi),qq(nxi,nqi),park(ns),pari(ni),
     *         rads(nh),par(kpars,nh,ncs),sfti1(nl,ncn),sfti2(nl,ncn),
     *          sip(nl,41),sih(nl,41),sib(nl,41),
     *          alfa(nl,ncn),beta(nl,ncn),sfti3(nl,ncn),
     *          psf(nl,ncn),pst(nl,ncn)
      data pi/3.14159265359/,re/6371.02e5/
      pmc=pi/180.
      msum(1)=0
      do1n=2,nqi
        nm=n-1
        msum(n)=msum(nm)+ntsl(nm)
    1 continue
      do 12 n=1,nqi
        nsum=msum(n)*kparn
        nt=ntsl(n)
        j1=ncn-n
        j2=n+1
        ls=0
        if(n.eq.1)goto3
          do 2 k=2,n
            ls=ls+ntsl(n-1)*2
    2     continue
   3   continue
        l=ls+1
          i1=1
          i3=(nt+1)/2
    5   continue
        i2=i1+1
        h=park(l)
        r=re+h
        ro=re/r
        ros=ro*ro
        roq=ros*ro
        t=park(l+1)*pmc
        ct=cos(t)
        st=sin(t)
        sks=1.+3.*ct*ct
        sk=sqrt(sks)
        sts=st*st
        b=bdip(h,t)
        q=qq(i1,n)
        ll=nsum+kparn*(i1-1)+1
        cm=pari(ll+3)
        co2=pari(ll)
        cn2=pari(ll+1)
        co=pari(ll+2)
        tn=pari(ll+4)
        vnu=pari(ll+6)
        vnv=pari(ll+7)
        call conducn(b,cm,co2,cn2,co,tn,sp,sh)
        s1=0.
        s2=0.
        s31=0.
        s32=0.
        s33=0.
        s4=0.
        s5=0.
        s9=0.
	s10=0.
	s11=0.
        alk2=st*sk/(r*r)
	bek2=1/(r*st)
        do7i=i2,i3
          l=l+2
          h=park(l)
          rp=re+h
          rop=re/rp
          rosp=rop*rop
          roqp=rosp*rop
          t=park(l+1)*pmc
          ctp=cos(t)
          stp=sin(t)
          sksp=1.+3.*ctp*ctp
          skp=sqrt(sksp)
          stsp=stp*stp
          bp=bdip(h,t)
          qp=qq(i,n)
          ll=nsum+kparn*(i-1)+1
          cmp=pari(ll+3)
          co2p=pari(ll)
          cn2p=pari(ll+1)
          cop=pari(ll+2)
          tnp=pari(ll+4)
          vnup=pari(ll+6)
          vnvp=pari(ll+7)
          call conducn(bp,cmp,co2p,cn2p,cop,tnp,spp,shp)
          dq=(qp-q)*.5
          s1=s1+(sp*r**3*st/(ro*sk)+spp*rp**3*stp/(rop*skp))*dq
          s2=s2+(sp*r*r/(roq*st*sks)+spp*rp*rp/(roqp*stp*sksp))*dq
          s31=s31+(sh*r**3/(ros*st*sks)+shp*rp**3/(rosp*stp*sksp))*dq
      	  s32=s32+(sh*r*r*st/(ros*sk)+shp*rp*rp*stp/(rosp*skp))*dq
          s33=s33+(sh*r*r/(roq*st*sk)+shp*rp*rp/(roqp*stp*skp))*dq*.5
	  s9=s9+(sp/(roq*sk)+spp/(roqp*skp))*dq*re
	  s10=s10+(sh/(roq*sk)+shp/(roqp*skp))*dq*re
          if(nvar.eq.2)goto6
            s4=s4+((sp*vnu-sh*vnv)*b*r*r/(roq*st*sks)+(spp*vnup-
     *      shp*vnvp)*bp*rp*rp/(roqp*stp*sksp))*dq
            s5=s5+((sh*vnu+sp*vnv)*b*r*r*st/(ros*sk)+
     *      (shp*vnup+spp*vnvp)*bp*rp*rp*stp/(rosp*skp))*dq
	      s11=s11+((sp*vnu-sh*vnv)*b/(roq*sk)+
     *	              (spp*vnup-shp*vnvp)*bp/(roqp*skp))*dq*re
    6     continue
          ct=ctp
          st=stp
          sks=sksp
          sk=skp
          sts=stsp
          r=rp
          ro=rop
          ros=rosp
          roq=roqp
          b=bp
          q=qp
          vnu=vnup
          vnv=vnvp
	  sp=spp
	  sh=shp
    7   continue
        alk1=st*sk/(r*r)
	bek1=1/(r*st)
        if(i1.ne.1)goto8
          alfa(kk,j1)=s1*alk1
          beta(kk,j1)=s2*bek1
          sfti1(kk,j1)=s31
          sfti2(kk,j1)=s32
          sfti3(kk,j1)=s33
          psf(kk,j1)=s4
          pst(kk,j1)=s5
	  sip(kk,j1)=s9
	  sih(kk,j1)=s10
	  sib(kk,j1)=s11
          goto9
    8   continue
          alfa(kk,j2)=s1*alk2
          beta(kk,j2)=s2*bek2
          sfti1(kk,j2)=s31
          sfti2(kk,j2)=s32
          sfti3(kk,j2)=s33
          psf(kk,j2)=s4
          pst(kk,j2)=s5
	  sip(kk,j2)=s9
	  sih(kk,j2)=s10
	  sib(kk,j2)=s11
    9   continue
        if(i3.eq.nt)goto12
          i1=(nt+2)/2
          if(i1.eq.i3)goto10
            l=l+2
   10     continue
          i3=nt
          goto5
   12 continue
      i=(ncn+1)/2
      pst(kk,i)=0.
      pst(kk,1)=0.
      pst(kk,ncn)=0.
      alfa(kk,i)=0.
      alfa(kk,1)=0.
      alfa(kk,ncn)=0.
      sfti1(kk,i)=0.
	  sfti2(kk,1)=0.
	  sfti2(kk,ncn)=0.
	  sfti2(kk,i)=0.
	  sfti3(kk,1)=0.
	  sfti3(kk,ncn)=0.
	  sfti3(kk,i)=0.
	  sip(kk,1)=0.
	  sih(kk,1)=0.
	  sib(kk,1)=0.
	  sip(kk,ncn)=0.
	  sih(kk,ncn)=0.
	  sib(kk,ncn)=0.
      return
      end
!------------------------------------------------------------------------------
      subroutine conducn(bb,cm,co2,cn2,co,tn,sp,sh)
      data e/1.60219e-20/,oe/1.76e7/,oi/3.09e2/,
     *ci1/4.23e-10/,ci2/4.28e-10/,ci3/2.58e-10/,
     *ce1/1.82e-10/,ce11/3.6e-2/,ce2/2.33e-11/,ce21/1.21e-4/,
     *ce3/2.8e-10/
      ee=e/bb*cm
      ome=oe*bb
      omi=oi*bb
      b=sqrt(tn)
      fi=ci1*co2+ci2*cn2+ci3*co
      fe=ce1*(1.+ce11* b)*b*co2+ce2*(1.-ce21*tn)*tn*cn2+ce3*b*co
      bu=fi/omi
      cu=fe/ome
      a=1./(1.+bu*bu)
      b=1./(1.+cu*cu)
      sp=ee*(bu*a+cu*b)
      sh=ee*(a-b)
      return
      end
!------------------------------------------------------------------------------
      subroutine jet_pB2(pmc,j1,nl,ddolgs,na,alt,ncn,dtet,nvar,
     *           pglo,kpar,nh,ncs,ids,pot0,curh,nzapjet,lzapjet)
         
      dimension alt(nh),pglo(kpar,nh,ncs,ids),pot0(na,nl,ncn),
     *          curh(nl,ncn,na),nzapjet(5),lzapjet(5)
      data re/6371.02e5/
      open(16,file='jet_p_B2',access='direct',recl=lzapjet(4))
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
               cm=pglo(6,k,j,i)
               co2=pglo(1,k,j,i)
               cn2=pglo(2,k,j,i)
               co=pglo(3,k,j,i)
               tn=pglo(7,k,j,i)
                 call conducn(b,cm,co2,cn2,co,tn,sp,sh)
                 s=0.
                 if(nvar.ne.2)then
                   si=2.*ct/sk
                   ci=st/sk
                   vnu=pglo(11,k,j,i)*si-pglo(10,k,j,i)*ci
                   vnv=pglo(12,k,j,i)
                   s=s+(sp*vnu-sh*vnv)*b
                 end if
                 efv=-(pot0(k,ip,j)-pot0(k,im,j))/(2.*ddolgs*pmc*r*st)
                 if(j.ne.j1)then
                    efu=-(pot0(k,i,j+1)-pot0(k,i,j-1))*sk/
     *			     (4.*dtet*pmc*r*ct)
                 else
                   efu=(pot0(kp,i,j)-pot0(km,i,j))/(alt(kp)-alt(km))
                 end if
                 curh(i,j,k)=efv*sp+efu*sh+s
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
      subroutine jetB2_80(pmc,j1,j2,r,alt,nh,nl,ncn,dtet,pot0,na,
     *  ddolgs,pef,nc,u,cur80,sip,sih,sib,nzapjet,lzapjet)
      dimension alt(nh),pot0(na,nl,ncn),pef(nl,nc),u(20),cur80(nl,nc),
     *  sip(nl,nc),sih(nl,nc),sib(nl,nc),nzapjet(5),lzapjet(5)
      data re/6371.02e5/
      open(14,file='jet_l_B2_80',access='direct',recl=lzapjet(2))
	do i=1,nl
	      ip=i+1
	      im=i-1
	      if(i.eq.1)im=nl
	      if(i.eq.nl)ip=1
	   do j=2,j1
	      js=nc-j+1
	      tet=(j-1)*dtet*pmc
	      tets=(js-5)*dtet*pmc
	      st=sin(tet)
	      sts=sin(tets)
	      ct=cos(tet)
	      cts=cos(tets)
	      sk=sqrt(1.+3.*ct*ct)
	      sks=sqrt(1.+3.*cts*cts)
	      efv=-(pot0(na,ip,j)-pot0(na,im,j))/(2.*ddolgs*pmc*r*st)
	      efvs=-(pot0(na,ip,js-4)-pot0(na,im,js-4))/
     *	                (2.*ddolgs*pmc*r*sts)
	      if(j.ne.j1)then 
	         efu=-(pot0(na,i,j+1)-pot0(na,i,j-1))*sk/
     *		      (4.*dtet*pmc*r*ct)
	         efus=-(pot0(na,i,js-3)-pot0(na,i,js-5))*sks/
     *      	      (4.*dtet*pmc*r*cts)
	      else 
                       efu=-(pef(i,j+1)-pot0(na,i,j-1))*re/
     *    	   (r*r*(u(19)-u(17)))
                       efus=efu
                    end if
                    cur80(i,j)=efv*sip(i,j)+efu*sih(i,j)+sib(i,j)
                    cur80(i,js)=efvs*sip(i,js)+
     *	        efus*sih(i,js)+sib(i,js)
         end do
         cur80(i,21)=0.
         efu=-(pef(i,21)-pef(i,19))*re/(r*r*(u(20)-u(18)))
         efv=-(pef(ip,20)-pef(im,20))/(2.*ddolgs*pmc*r)
         cur80(i,20)=efv*sip(i,20)+efu*sih(i,20)+sib(i,20)
         cur80(i,22)=efv*sip(i,22)+efu*sih(i,22)+sib(i,22)
      end do
      jstart=2
      jfinish=nc-jstart+1
      do i=1,nl
         phi=(i-1)*ddolgs
         do j=jstart,jfinish
            if(j.eq.j2)then
               tet=90.
            else
               if(j.lt.j2)then
                  st=sqrt(u(j-1)*(re+alt(1))/re)
               else
                  js=nc-j+1
                  st=sqrt(u(js-1)*(re+alt(1))/re)
               end if
               tet=asin(st)/pmc
               if(j.gt.j2)tet=180.-tet
            end if
                nzapjet(2)=nzapjet(2)+1
                write(14,rec=nzapjet(2))phi,90.-tet,cur80(i,j)*1.e6,
     *          sip(i,j)*1.e9,sih(i,j)*1.e9
         end do
      end do
      close(14)
      return
      end
!------------------------------------------------------------------------------
      subroutine jetB2_h(pmc,j1,r,alt,nl,ncn,dtet,nvar,pot0,na,
     *           ddolgs,pglo,kpar,nh,ncs,ids,cur,nzapjet,lzapjet)
      dimension alt(nh),pot0(na,nl,ncn),cur(nl,ncn),
     *          pglo(kpar,nh,ncs,ids),nzapjet(5),lzapjet(5)
      open(15,file='jet_l_B2_h',access='direct',recl=lzapjet(3))
      ki=na
      kip=ki+1
      kim=ki-1
      if(ki.eq.1)kim=ki
      if(ki.eq.na)kip=ki
      do i=1,nl
         ip=i+1
         im=i-1
         if(i.eq.1)im=nl
         if(i.eq.nl)ip=1
         do j=2,j1
               tet=(j-1)*dtet*pmc
               h=alt(1)
               b=bdip(h,tet)
               cm=pglo(6,1,j,i)
               co2=pglo(1,1,j,i)
               cn2=pglo(2,1,j,i)
               co=pglo(3,1,j,i)
               tn=pglo(7,1,j,i)
               ct=cos(tet)
               st=sin(tet)
               sk=sqrt(1.+3.*ct*ct)
            if(nvar.ne.2)then
               si=2.*ct/sk
               ci=st/sk
               vnu=pglo(11,1,j,i)*si-pglo(10,1,j,i)*ci
               vnv=pglo(12,1,j,i)
            end if
               call conducn(b,cm,co2,cn2,co,tn,sp,sh)
               s=0.
               spi=0.
               shi=0.
               if(j.ne.j1)then
	    js=ncn-j+1
                  tets=(js-1)*dtet*pmc
                  bs=bdip(h,tets)
                  cms=pglo(6,1,js,i)
                  co2s=pglo(1,1,js,i)
                  cn2s=pglo(2,1,js,i)
                  coss=pglo(3,1,js,i)
                  tns=pglo(7,1,js,i)
                  cts=cos(tets)
                  sts=sin(tets)
                  sks=sqrt(1.+3.*cts*cts)
            if(nvar.ne.2)then
                  sis=2.*cts/sks
                  cis=sts/sks
                  vnus=pglo(11,1,js,i)*sis-pglo(10,1,js,i)*cis
                  vnvs=pglo(12,1,js,i)
            end if
                  call conducn(bs,cms,co2s,cn2s,coss,tns,sps,shs)
                  ss=0.
                  spis=0.
                  shis=0.
               end if
               do k=1,na-1
                  kp=k+1
                  hp=alt(kp)
                  bp=bdip(hp,tet)
                  cmp=pglo(6,kp,j,i)
                  co2p=pglo(1,kp,j,i)
                  cn2p=pglo(2,kp,j,i)
                  cop=pglo(3,kp,j,i)
                  tnp=pglo(7,kp,j,i)
            if(nvar.ne.2)then
                  vnup=pglo(11,kp,j,i)*si-pglo(10,kp,j,i)*ci
                  vnvp=pglo(12,kp,j,i)
            end if
                  call conducn(bp,cmp,co2p,cn2p,cop,tnp,spp,shp)
                  spi=spi+(sp+spp)*.5*(hp-h)
                  shi=shi+(sh+shp)*.5*(hp-h)
            if(nvar.ne.2)then
                  s=s+(sp*vnu*b+spp*vnup*bp)*.5*(hp-h)-
     *                   (sh*vnv*b+shp*vnvp*bp)*.5*(hp-h)
            end if
                  if(j.ne.j1)then
                     bsp=bdip(hp,tets)
                     cmsp=pglo(6,kp,js,i)
                     co2sp=pglo(1,kp,js,i)
                     cn2sp=pglo(2,kp,js,i)
                     cossp=pglo(3,kp,js,i)
                     tnsp=pglo(7,kp,js,i)
            if(nvar.ne.2)then
                     vnusp=pglo(11,kp,js,i)*sis-pglo(10,kp,js,i)*cis
                     vnvsp=pglo(12,kp,js,i)
            end if
                     call conducn(bsp,cmsp,co2sp,cn2sp,cossp,
     *			   tnsp,spsp,shsp)
                     spis=spis+(sps+spsp)*.5*(hp-h)
                     shis=shis+(shs+shsp)*.5*(hp-h)
            if(nvar.ne.2)then
                     ss=ss+(sps*vnus*bs+spsp*vnusp*bsp)*.5*(hp-h)-
     *                   (shs*vnvs*bs+shsp*vnvsp*bsp)*.5*(hp-h)
            end if
                     bs=bsp
                     cms=cmsp
                     co2s=co2sp
                     cn2s=cn2sp
                     coss=cossp
                     tns=tnsp
            if(nvar.ne.2)then
                     vnus=vnusp
                     vnvs=vnvsp
            end if
                     sps=spsp
                     shs=shsp                 
                  end if
                  h=hp
                  b=bp
                  cm=cmp
                  co2=co2p
                  cn2=cn2p
                  co=cop
                  tn=tnp
            if(nvar.ne.2)then
                  vnu=vnup
                  vnv=vnvp
            end if
                  sp=spp
                  sh=shp
               end do        
               efv=-(pot0(ki,ip,j)-pot0(ki,im,j))/(2.*ddolgs*pmc*r*st)
               if(j.ne.j1)then
	         efvs=-(pot0(ki,ip,js)-pot0(ki,im,js))/
     *	                (2.*ddolgs*pmc*r*sts)
	         efu=-(pot0(ki,i,j+1)-pot0(ki,i,j-1))*sk/
     *		      (4.*dtet*pmc*r*ct)
	         efus=-(pot0(ki,i,js+1)-pot0(ki,i,js-1))*sks/
     *                 	      (4.*dtet*pmc*r*cts)
               else 
                  efu=(pot0(kip,i,j)-pot0(kim,i,j))/(alt(kip)-alt(kim))
               end if
               cur(i,j)=efv*spi+efu*shi+s
               if(j.ne.j1)cur(i,js)=efvs*spis+efus*shis+ss
          end do
      end do
      jstart=2
      jfinish=ncn-jstart+1	 
      do i=1,nl
	   phi=(i-1)*ddolgs
	   do j=jstart,jfinish
	     tet=(j-1)*dtet
                   nzapjet(3)=nzapjet(3)+1
                   write(15,rec=nzapjet(3))phi,90.-tet,cur(i,j)*1.e6
	   end do
	end do
      close(15)
      return
      end
!------------------------------------------------------------------------------
      subroutine jetB2(pmc,j1,r,nl,ncn,dtet,pot0,na,ddolgs,pef,
     *           nc,u,cur,sip,sih,sib,nzapjet,lzapjet)
      dimension pot0(na,nl,ncn),pef(nl,nc),u(20),cur(nl,ncn),
     *          sip(nl,nc),sih(nl,nc),sib(nl,nc),nzapjet(5),lzapjet(5)
      data re/6371.02e5/
      open(13,file='jet_l_B2',access='direct',recl=lzapjet(1))
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
	      ct=cos(tet)
	      cts=cos(tets)
	      sk=sqrt(1.+3.*ct*ct)
	      sks=sqrt(1.+3.*cts*cts)
	      efv=-(pot0(na,ip,j)-pot0(na,im,j))/(2.*ddolgs*pmc*r*st)
	      if(j.ne.j1)then
                       efvs=-(pot0(na,ip,js)-pot0(na,im,js))/
     *	                (2.*ddolgs*pmc*r*sts)
	         efu=-(pot0(na,i,j+1)-pot0(na,i,j-1))*sk/
     *		      (4.*dtet*pmc*r*ct)
	         efus=-(pot0(na,i,js+1)-pot0(na,i,js-1))*sks/
     *      	      (4.*dtet*pmc*r*cts)
	      else 
               efu=-(pef(i,j+1)-pot0(na,i,j-1))*re/
     *        		 (r*r*(u(19)-u(17)))
            end if
            cur(i,j)=efv*sip(i,j)+efu*sih(i,j)+sib(i,j)
            if(j.ne.j1)cur(i,js)=efvs*sip(i,js+4)+
     *                 efus*sih(i,js+4)+sib(i,js+4)
            end do
	end do 
	jstart=2
	jfinish=ncn-jstart+1	 
        do i=1,nl
	   phi=(i-1)*ddolgs
	   do j=jstart,jfinish
             nzapjet(1)=nzapjet(1)+1
	     tet=(j-1)*dtet
	     if(j.lt.j1)then
                 write(13,rec=nzapjet(1))phi,90.-tet,cur(i,j)*1.e6,
     *           sip(i,j)*1.e9,sih(i,j)*1.e9
	     else if(j.eq.j1)then
                 write(13,rec=nzapjet(1))phi,90.-tet,cur(i,j)*1.e6,
     *          (sip(i,j)+sip(i,j+4))*1.e9,	        
     *          (sih(i,j)+sih(i,j+4))*1.e9
	     else
              write(13,rec=nzapjet(1))phi,90.-tet,cur(i,j)*1.e6,
     *         sip(i,j+4)*1.e9,	        
     *         sih(i,j+4)*1.e9
	     end if
          end do
      end do
      close(13)
      return 
      end
!------------------------------------------------------------------------------
      subroutine degapn(nvar,nc,nl,dtet,dfi,sfti1,sfti2,sfti3,
     *            psf,pst,gamma,delta,psi,u,nui,rads,nh,na)
      dimension rads(nh),gamma(nl,nc),delta(nl,nc),psi(nl,nc),
     *          pst(nl,nc),psf(nl,nc),sfti1(nl,nc),u(nui),sfti2(nl,nc),
     *          sfti3(nl,nc)
      data pi/3.14159265359/,re/6371.02e5/
      df=(dfi+dfi)/180.*pi
      dt=(dtet+dtet)/180.*pi
      df=1./df
      dt=1./dt
      r=re+rads(na)
      s1=(re+re)/r
c      s1=re/r
      k=nc-1
      k1=(nc+1)/2
      k2=k1-3
      k3=k1+3
      do5i=2,k
        if(i.eq.k1)goto5
          ip=i+1
          im=i-1
          if(i.le.k2)tet=im*dtet/180.*pi
          if(i.ge.k3)tet=(i-5)*dtet/180.*pi
          if(i.gt.k2.and.i.lt.k3)goto1
	      st=sin(tet)
	      ct=cos(tet)
	      sks=sqrt(1.+3.*ct*ct)
            s=s1*st*ct
c            ct=cos(tet)
c            s=s1*(1.+3.*ct*ct)*sin(tet)*.5/ct
            goto2
    1     continue
            if(i.lt.k1)l=i-1
            if(i.lt.k1)goto2
              l=nc-i
    2       continue
            do4j=1,nl
              jp=j+1
              jm=j-1
              if(j.eq.1)jm=nl
              if(j.eq.nl)jp=1
              gamma(j,i)=(sfti1(jp,i)-sfti1(jm,i))*df*st*sks/(r*r)
              psi(j,i)=0.
              if(nvar.ne.2)psi(j,i)=(psf(jp,i)-psf(jm,i))*df
              if(i.gt.k2.and.i.lt.k3) goto3
              delta(j,i)=(sfti2(j,ip)-sfti2(j,im))*dt/(s*r*st)+
     *                    sfti3(j,i)/(r*st) 
              if(nvar.ne.2) psi(j,i)=psi(j,i)-(pst(j,ip)-pst(j,im))*
     *                      dt/s
              goto4
    3         continue
              du=u(l+1)-u(l-1)
              if(i.gt.k1)du=-du
              delta(j,i)=(sfti2(j,ip)-sfti2(j,im))/(du*r*st)+
     *                    sfti3(j,i)/(r*st)
              if(nvar.ne.2) psi(j,i)=psi(j,i)-(pst(j,ip)-pst(j,im))/du
    4       continue
    5 continue
c     do6i=1,nl
c       do6j=5,k2
c         l=nc-j+1
c         gamma(i,j)=gamma(i,j)+gamma(i,l)
c         gamma(i,l)=gamma(i,j)
c         delta(i,j)=delta(i,j)+delta(i,l)
c         delta(i,l)=delta(i,j)
c         psi(i,j)=psi(i,j)+psi(i,l)
c         psi(i,l)=psi(i,j)
c   6 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine sftpoln(ll,rads,nh,na,par,kpars,ncs,ncn,sfti,nl)
      dimension rads(nh),par(kpars,nh,ncs),sfti(nl,ncn)
      sfti(ll,1)=0.
      sfti(ll,ncn)=0.
      return
      end
!------------------------------------------------------------------------------
      subroutine libmanmf(nc,nl,na,alt,nh,dtets,dfis,alfa,
     *           beta,gamma,delta,psi,pef0,om,tet1,tet2,tet3,
     *           ae,be,ge,pef,u,nui,nvg)
      dimension alt(nh),alfa(nl,nc),beta(nl,nc),gamma(nl,nc),ae(nl),
     *          delta(nl,nc),psi(nl,nc),pef0(nl,nc),pef(nl,nc),be(nl),
     *          ge(nl),u(nui)
      data pi/3.14159265359/,re/6371.02e5/,pmi0/14.9/
      df=dfis/180.*pi
      dt=dtets/180.*pi
      dt1=dt*.5
      df=1./df
      dt=1./dt
      k=nc-1
      dfs=df*df*.5
      dts=dt*dt*.5
      df=df*.5
      dt=dt*.5
      m=(nc+1)/2
      mm=m-3
      mp=m+3
c      s1=(re+alt(na))/re*2.
      s1=(re+alt(na))/(re*2.)
      s2=s1*s1
      do4i=1,nl
        ip=i+1
        im=i-1
        if(i.eq.1)im=nl
        if(i.eq.nl)ip=1
        do3j=2,k
          if(j.eq.m)then
            call libme(i,ip,im,j,u,nui,nl,ae,be,ge,nc,pef0,pef,om,df)
c           pef(i,j)=0.
c           pef(i,j)=pef(i,j-1)
          else
            p=0.
            f=0.
            if(j.gt.mm.and.j.lt.mp)then
              l=nc-j+1
              n=j-1
              if(j.lt.m)then
                call libmu(i,ip,im,j,l,n,nl,nc,alfa,beta,gamma,delta,
     *          psi,u,nui,df,dfs,pef0,pef,om,p,f)
              else
                pef(i,j)=pef(i,l)
              end if
              goto3
            else
              if(j.le.mm)then
                tetn=(j-1)*dtets/180.*pi
              else
                tetn=(j-5)*dtets/180.*pi
              end if
              st=sin(tetn)
              sts=st*st
              pmi=(re+alt(na))/(re*sts)
              if(pmi.ge.pmi0)goto1
                l=nc-j+1
                if(j.gt.m)then
                  pef(i,j)=pef(i,l)
                  goto3
                else
c                  if(l.le.mm)tets=(l-1)*dtets/180.*pi
c                  if(l.ge.mp)tets=(l-5)*dtets/180.*pi
                  tets=(l-5)*dtets/180.*pi
                  call libmsn(i,ip,im,l,alfa,beta,gamma,delta,psi,dt,
     *            dts,df,dfs,pef0,nl,nc,p,f,s1,s2,dt1,tets)
c    *            dts,df,dfs,pef,nl,nc,p,f,s1,s2,dt1,tets)
                end if
    1         continue
              if(nvg.eq.0)goto2
                if(j.eq.4.or.j.eq.38)goto3
    2         continue
              call libmsn(i,ip,im,j,alfa,beta,gamma,delta,psi,dt,
     *        dts,df,dfs,pef0,nl,nc,p,f,s1,s2,dt1,tetn)
c    *        dts,df,dfs,pef,nl,nc,p,f,s1,s2,dt1,tetn)
              pef(i,j)=om/p*f+(1.-om)*pef0(i,j)
            end if
          end if
    3   continue
    4 continue
      s1=0.
      s2=0.
      do5i=1,nl
        s1=s1+pef(i,2)
        s2=s2+pef(i,k)
    5 continue
      s1=s1/nl
      s2=s2/nl
      do6i=1,nl
        pef(i,1)=s1
        pef(i,nc)=s2
    6 continue
      return
      end
!------------------------------------------------------------------------------
      subroutine libme(i,ip,im,j,u,nui,nl,ae,be,ge,nc,pef0,pef,om,df)
      dimension u(nui),ae(nl),be(nl),ge(nl),pef0(nl,nc),pef(nl,nc)
      jm=j-1
      a=1./(u(nui)-u(nui-1))
      b=-be(i)*a
      e=ae(i)
      a=abs(e)
      g=a+e
      e=e-a
c     pef(i,j)=1./(om*(b+(a+a)*df)+1.)*((b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))*om+pef0(i,j))
c     pef(i,j)=om/(om*(b+(a+a)*df)+4.)*((b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))*2.+(4.-om*(b+(a+a)*df))/om*pef0(i,j))
c     pef(i,j)=(.5*om/(b+(a+a)*df)*(b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))+(1.-om*.25)*pef0(i,j))/(1.+.25*om)
c     pef(i,j)=om/(b+(a+a)*df)*(b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))+(1.-om)*pef0(i,j)
      pef(i,j)=om/(b+(a+a)*df)*(b*pef0(i,jm)+(g*pef0(ip,j)-
     *e*pef0(im,j))*df-ge(i))+(1.-om)*pef0(i,j)
      return
      end
!------------------------------------------------------------------------------
      subroutine libmsn(i,ip,im,j,alfa,beta,gamma,delta,psi,
     *dt,dts,df,dfs,pef0,nl,nc,p,f,s1,s2,dt1,tet)
      dimension alfa(nl,nc),beta(nl,nc),gamma(nl,nc),
     *delta(nl,nc),psi(nl,nc),pef0(nl,nc)
      jp=j+1
      jm=j-1
      e=alfa(i,j)
      g=(e+alfa(i,jp))*dts
      s3=sin(tet)*cos(tet)
c      ct=cos(tet)
c      s3=(1.+3.*ct*ct)*sin(tet)/ct
      tets=tet+dt1
      s4=sin(tets)*cos(tets)
c      ct=cos(tets)
c      s4=(1.+3.*ct*ct)*sin(tets)/ct
      s5=s2/s3
      g=g*s5/s4
      h=gamma(i,j)*dt*s1/s3
c      h=gamma(i,j)*dt*s1/s4
      a=abs(h)
      b=g+h+a
      p=p+b
      f=f+b*pef0(i,jp)
      g=(e+alfa(i,jm))*dts
      tets=tet-dt1
      s4=sin(tets)*cos(tets)
c      ct=cos(tets)
c      s4=(1.+3.*ct*ct)*sin(tets)/ct
      g=g*s5/s4
      b=g-h+a
      p=p+b
      f=f+b*pef0(i,jm)
      e=beta(i,j)
      h=delta(i,j)*df
      g=(e+beta(ip,j))*dfs
      a=abs(h)
      b=g-h+a
      p=p+b
      f=f+b*pef0(ip,j)
      g=(e+beta(im,j))*dfs
      b=g+h+a
      p=p+b
      f=f+b*pef0(im,j)
      f=f-psi(i,j)
      return
      end
!------------------------------------------------------------------------------
      subroutine libmu(i,ip,im,j,l,n,nl,nc,alfa,beta,gamma,delta,psi,
     *u,nui,df,dfs,pef0,pef,om,p,f)
      dimension alfa(nl,nc),beta(nl,nc),gamma(nl,nc),u(nui),
     *delta(nl,nc),psi(nl,nc),pef0(nl,nc),pef(nl,nc)
      jp=j+1
      jm=j-1
      lp=l-1
      lm=l+1
      du=u(n+1)-u(n)
      du=1./du
      dus=du*du*.5
      du=du*.5
      e=alfa(i,l)
      g=(e+alfa(i,lp))*dus
      h=gamma(i,l)*du
      a=abs(h)
      b=g+h+a
      p=p+b
c     f=f+b*pef(i,lp)
      f=f+b*pef0(i,lp)
      g=(e+alfa(i,lm))*dus
      b=g-h+a
      p=p+b
c     f=f+b*pef(i,lm)
      f=f+b*pef0(i,lm)
      e=beta(i,l)
      h=delta(i,l)*df
      g=(e+beta(ip,l))*dfs
      a=abs(h)
      b=g-h+a
      p=p+b
c     f=f+b*pef(ip,l)
      f=f+b*pef0(ip,l)
      g=(e+beta(im,l))*dfs
      b=g+h+a
      p=p+b
c     f=f+b*pef(im,l)-psi(i,l)
      f=f+b*pef0(im,l)-psi(i,l)
      e=alfa(i,j)
      h=gamma(i,j)*du
      g=(e+alfa(i,jp))*dus
      a=abs(h)
      b=g+h+a
      p=p+b
c     f=f+b*pef(i,jp)
      f=f+b*pef0(i,jp)
      g=(e+alfa(i,jm))*dus
      b=g-h+a
      p=p+b
c     f=f+b*pef(i,jm)
      f=f+b*pef0(i,jm)
      e=beta(i,j)
      h=delta(i,j)*df
      g=(e+beta(ip,j))*dfs
      a=abs(h)
      b=g-h+a
      p=p+b
c     f=f+b*pef(ip,j)
      f=f+b*pef0(ip,j)
      g=(e+beta(im,j))*dfs
      b=g+h+a
      p=p+b
c     f=f+b*pef(im,j)
      f=f+b*pef0(im,j)
      f=f-psi(i,j)
      pef(i,j)=om/p*f+(1.-om)*pef0(i,j)
      return
      end
!------------------------------------------------------------------------------