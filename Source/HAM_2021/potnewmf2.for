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
