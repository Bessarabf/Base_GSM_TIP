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