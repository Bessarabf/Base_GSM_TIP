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