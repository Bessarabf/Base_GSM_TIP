      subroutine jetB1_h(pmc,j1,r,alt,nl,ncn,dtet,nvar,pot0,na,
     *           ddolgs,pglo,kpar,nh,ncs,ids,cur,nzapjet,lzapjet)
      dimension alt(nh),pot0(na,nl,ncn),cur(nl,ncn),
     *          pglo(kpar,nh,ncs,ids),nzapjet(5),lzapjet(5)
      open(15,file='jet_l_B1_h',access='direct',recl=lzapjet(3))
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