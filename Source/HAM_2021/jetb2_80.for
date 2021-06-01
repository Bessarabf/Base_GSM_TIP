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