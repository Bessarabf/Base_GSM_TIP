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
