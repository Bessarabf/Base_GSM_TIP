      subroutine iseti(nh,na,nui,u,nxi,nqi,q,rads,ns,park,ntsl)
      dimension rads(nh),park(*),u(*),q(32,19),ntsl(*)
      double precision pi,a,b,c,d,f,g,p,re,e
      rmin=80.
      hm=rmin*1.e5
      nxi=na+na
      dh=3.e5
      gamma=1.1
      rads(1)=hm
      do1i=2,nh
        im=i-1
        rads(i)=rads(im)+dh*gamma**(im-1)
    1 continue
      ns=0
      e=rads(1)
      nl=0
      b=175.d0
      c=90.d0
      re=6.37102d8 ! double presision 17.10.2013
      pi=3.14159265359d0
      a=pi/180.
      d=-5.d0
      htm=rads(na)
      il=0
    2 continue
      if(b.lt.c)go to 11
    3   continue
        nl=nl+1
        li=ns+1
        park(li)=e
        f=dsin(b*a)
        f=f*f/(re+htm)
        if(il.ne.0)f=1.d0/re*(u(nl-1)*2.-u(nl-2))
        p=(re+e)*f
        park(li+1)=pi-datan(dsqrt(p/(1.d0-p)))
        g=e+dh
        h=1./f-re
        o=amin1(h,htm)
        li=li+2
        i=1
    4   continue
        if(g.ge.o)go to 5
          i=i+1
          park(li)=g
          p=(re+g)*f
          park(li+1)=pi-datan(dsqrt(p/(1.d0-p)))
          g=g+dh*gamma**(i-1)
          li=li+2
          go to 4
    5   continue
        j=i+1
        if(o.ne.h)go to 7
          park(li)=h
          park(li+1)=pi*.5
          nx=j+j-1
          do 6 k=1,i
            m=ns+2*(j+k)-1
            l=ns+2*(j-k)-1
            park(m)=park(l)
            park(m+1)=pi-park(l+1)
    6     continue
          go to 9
    7   continue
          park(li)=o
          g=(re+o)*f
          h=pi-datan(dsqrt(g/(1.d0-g)))
          park(li+1)=h
c         li=li+2
c         park(li)=o
c         park(li+1)=pi-h
          nx=j+j
          do 8 k=1,j
            l=ns+2*(j-k)+1
            m=ns+2*(j+k)-1
            park(m)=park(l)
            park(m+1)=pi-park(l+1)
    8     continue
    9   continue
        li=ns+1
        f= sin(park(ns+2))
        u(nl)=re/(re+park(li))*f*f
        do 10 i=1,nx
          o=park(li)
          f=park(li+1)
          h=re/(re+o)
          q(i,nl)=h*h*dcos(f)
          li=li+2
   10   continue
        i=nx/2
        j=(nx+1)/2
        if(i.ne.j) q(j,nl)=0.
        ntsl(nl)=nx
        ns=ns+nx*2
        if(il.eq.1)goto12
          b=b+d
          go to 2
   11 continue
      il=1
      goto3
   12   continue
        do 13 i=2,ns,2
          park(i)=park(i)/a
   13   continue
        nqi=nl
        nui=nl+1
        u(nui)=re/(re+hm)
        return
      end
