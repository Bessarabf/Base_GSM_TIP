      subroutine field(et,el,vet,vel,its,
     *           dtet,ddol,potef,ntr,nl2,ids,i5,j)
      dimension
     *          potef(ntr,ids,nl2)
      data b0/0.3/,
     *     pi/3.14159/,re/6.371e 8/
      cr=180./pi
      itd=its-1
      dtet5=dtet*.5
      tet=dtet5*(i5-1)
      ct=cos(tet)
      st=sin(tet)
      bb=b0*sqrt(1+3.*ct*ct)
      dd2=2.*ddol
      dt2=2.*dtet5
      jp=j+1
      jm=j-1
      if(j.eq.1) jm=ids
      if(j.eq.ids) jp=1
      def=(potef(ntr,jp,i5)-potef(ntr,jm,i5))/dd2
      if(i5.eq.4) then
       det=(potef(ntr,j,i5)-potef(ntr,j,i5-1))/dtet5
      else
       det=(potef(ntr,j,i5+1)-potef(ntr,j,i5))/dtet5
      end if
cc    det=(potef(ntr,j,i5+1)-potef(ntr,j,i5-1))/dt2
      et=-det/re
      el=-def/st/re
      vel=el/bb
      vet=et/bb
      return
      end
