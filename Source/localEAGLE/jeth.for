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