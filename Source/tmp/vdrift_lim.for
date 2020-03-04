! 28/02/2020 limit on dreif velocity
      subroutine vdrift(i,j,idt,ntsl,nl,ddolg,dtet,fi,ht,tt,
     *           potef,vdu,vdv,ntr,nl2)

      dimension tt(*),ntsl(nl),vdu(*),vdv(*),
     *          potef(ntr,idt,nl2),ht(*)

      data re/6371.02e5/,pi/3.14159265359/
      im=i-1
      if(i.eq.1)im=idt
      ip=i+1
      if(i.eq.idt)ip=1
      dt=dtet/90.*pi
      dd=ddolg/90.*pi
      k2=(ntsl(j)+1)/2
      r=re+ht(1)
      tp=tt(1)
      k1=1
      jp=nl2-j
      jm=jp+1
      jd=jp-1
      st=sin(tp)
      ct=cos(tp)
    1 continue
      dff=(potef(16,ip,jp)-potef(16,im,jp))/dd
      dfu=(potef(16,i,jm)-potef(16,i,jd))/dt
      if(jp.eq.4.or.jp.eq.33)dfu=(potef(16,i,jp)-potef(16,i,jd))/dt*2.
      if(jp.eq.5.or.jp.eq.34)dfu=(potef(16,i,jm)-potef(16,i,jp))/dt*2.
! limit on dreif velocity     
      dfu=dfu/r
      if(abs(dfu).gt.1.e5) then
        print*,'WARNING: dreif velocity too large', dfu,i,j
        dfu=dfu/abs(dfu)*1.e5
      end if
      dfu=dfu*.5*st*st*st/ct
! lim end
      dff=dff*st*st/r
      do 2 k=k1,k2
        h=ht(k)
        t=tt(k)
        if(k.eq.ntsl(j))t=tt(1)
        b=bdip(h,t)
        st=sin(t)
        ct=cos(t)
        cev=1./(b*st**3)
        vdu(k)=dff*cev
        vdv(k)=-dfu*sqrt(1.+3.*ct*ct)*cev
    2 continue
      if(k1.ne.1)go to 3
         k1=k2+1
         k2=ntsl(j)
         tp=tt(1)
c        tp=tt(k2)
         jp=j+1
         jd=j
         jm=jp+1
         st=sin(tp)
         ct=-cos(tp)
c        ct= cos(tp)
         go to 1
    3 continue
c     print4,vdu(1),vdu(k2),vdv(1),vdv(k2)
c   4 format(' ',3x,'vu(1)=',g12.3,3x,'vu(n)=',g12.3,3x,
c    *'vv(1)=',g12.3,3x,'vv(n)=',g12.3)
      return
      end