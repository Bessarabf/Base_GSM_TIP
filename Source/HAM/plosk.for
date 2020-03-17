c
      subroutine plosk(pk,pk3,hm,cmi,ssh,se,cm1,ns,ip,imin,ntsl
     *   ,ntpl,nl,npar,its,nh,nl2)

      integer ntsl(nl)
      real hm(ip),cm1(ip,ns),pk(8,ntpl),
     *   pk3(2,ntpl),cmi(its,nh),ssh(its),se(ns)
	allocatable sm(:),cm(:)
      allocate (sm(nl2-2),cm(nl2-2))      
	ih=imin
      do 29 ihp=1,ip
      hm(ihp)=pk3(1,ih)*1.e-5
      i0=0
      i=0
      do 30 jt=1,nl
      n2=ntsl(jt)/2
      if(ih.gt.n2)goto 19
      i=i+2
      ih1=i0+ih
      sm(i)=pk3(2,ih1)
      cm(i)=cf(pk,ntpl,8,npar,ih1)
      i0=i0+ntsl(jt)
      ih1=i0+1-ih
      sm(i-1)=pk3(2,ih1)
      cm(i-1)=cf(pk,ntpl,8,npar,ih1)
  30  continue
  19  call turn1(sm,cm,i)
c      call spline(i,sm,cm,b,c,d)
      do 33 j=1,ns
      cm1(ihp,j)=val(i,se(j),sm,cm)
      if(ih.le.15.and.npar.eq.9) cm1(ihp,j)=alog10(10.**(cm1(ihp,j))+
     *        10.**val(19,se(j),ssh,cmi(1,ih+15)))
  33  continue
  29  ih=ih+1
      deallocate (sm,cm)
      return
      end
c
      function cf(pk,ntpl,npk,npar,j1)
      real pk(npk,ntpl)
      if(npar.ne.9)then
        cf=pk(npar,j1)
        if(npar.ge.4.and.npar.le.6) cf=cf*pk(npar-3,j1)
        if(npar.le.3)then
c       if(npar.le.3.or.npar.ge.7)then
           cf=alog10(cf)
cc         if(cf.lt.0.) cf=0.
           if(cf.lt.-3.) cf=-3.
        endif
      else
        cf=alog10(pk(1,j1)+pk(2,j1))
      endif
      return
      end
c
      subroutine turn1(h,c,n)
      real h(n),c(n)
      j=1
  3   if(h(j).le.h(j+1))goto 2
      call turn(h(j))
      call turn(c(j))
      j=j-2
      if(j.eq.-1)j=0
  2   j=j+1
      if(j.lt.n)goto 3
      return
      end
c
      subroutine turn(sm)
      real sm(2)
      a=sm(1)
      sm(1)=sm(2)
      sm(2)=a
      return
      end
cc
      real function val(n,u,x,y)
      real x(n),y(n)
      data i/1/,k/1/
      if(i.ge.n)i=1
      if(u.lt.x(i)) go to 10
      if(u.le.x(i+1)) go to 30
   10 i=1
      j=n+1
   20 k=(i+j)/2
      if( u.lt.x(k) ) j=k
      if(u.ge.x(k)) i=k
      if(j.gt.(i+1)) go to 20
   30 dx=(u-x(i))/(x(i+1)-x(i))
      val=y(i)+dx*(y(i+1)-y(i))
      return
      end


