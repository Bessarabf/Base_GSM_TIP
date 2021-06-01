c                S E T K A
c ver 02.03.2018 
      subroutine setka_kl(rmin,dh,gamma,ntr,nh,b1,c1,dtett,
     *                 rmaxt,par,nr,gkoor,nv,
     *                 rads,ntsl,nl,q,u,its,ids,dtets,ddolgs)
      dimension par(nr),rads(nh),
     *          q(nv,nl) ,u(nl),ntsl(nl),gkoor(2,its,ids)
      double precision pi,a,b,c,d,f,g,p,re,e
      data  pi/3.14159265359d0/
      b=b1
      c=c1
      re=6.37102d8
      a=pi/180.d0
      d=-dtett
      htm=rmaxt*1.e5-re
      hm=rmin*1.e5
      dh=dh*1.e5
      rads(1)=hm
c  . . . высоты "шара"
      do i=2,nh
        rads(i)= rads(i-1)+dh*gamma**(i-2)
      end do
      ns=0
      e=rads(ntr)
      nl=0
C    . . . сетка "трубки":
C    . . . цикл по коширотам от юж.полюса до экватора
    2 continue
      if(b.lt.c)go to 14
        nl=nl+1
        li=ns+1
        par(li)=e
        f=dsin(b*a)
        f=f*f/(re+e)
        par(li+1)=b*a
        g=e+dh*gamma**(ntr-1)
        h=1.d0/f-re           ! . . . высота вершины линии
        o=amin1(h,htm)        ! . . . высота верхнего узла линии
        li=li+2
        i=1
c . . . Цикл по I с предусловием
    3   if(g.gt.o)go to  4
C   . . .  Восходящая половина линии
          i=i+1
          par(li)=g
          p=(re+g)*f
          par(li+1)=pi-datan(dsqrt(p/(1.d0-p)))
          g=g+dh*gamma**(i+ntr-2)
          li=li+2
          go to 3
    4   continue
c . . . Конец цикла по I
        j=i+1
        if(o.eq.h) then
c . . . Замкнутая линия:
          par(li)=h              ! . . . экваториальная точка
          par(li+1)=pi*.5
          nx=j+j-1
          do k=1,i
            m=ns+2*(j+k)-1
            l=ns+2*(j-k)-1
            par(m)=par(l)
            par(m+1)=pi-par(l+1)
! ...широты основания линий - одинаковы! 02.03.2018
            if(k.eq.i)par(m+1)=(180.-b)*a 
          end do
        else
c . . . Разомкнутая линия:
          par(li)=o
          g=(re+o)*f
          h=pi-datan(dsqrt(g/(1.d0-g)))
          par(li+1)=h
          li=li+2
          par(li)=o
          par(li+1)=pi-h
          nx=j+j
          do k=1,j
            l=ns+2*(j-k)+1
            m=ns+2*(j+k)-1
            par(m)=par(l)
            par(m+1)=pi-par(l+1)
! ...широты основания линий - одинаковы! 02.03.2018
            if(k.eq.j)par(m+1)=(180.-b)*a
          end do
        end if
        li=ns+1
        f= sin(par(ns+2))
        u(nl)=re/(re+par(li))*f*f
        do i=1,nx
          o=par(li)
          f=par(li+1)
          h=re/(re+o)
          q(i,nl)=h*h*dcos(f)
          li=li+2
        end do
        i=nx/2
        j=(nx+1)/2
        if(i.ne.j) q(j,nl)=0.
        ntsl(nl)=nx
        ns=ns+nx*2
        b=b+d
        go to 2
   14 continue
c . . . Преревод в градусы коширот
      do i=2,ns,2
        par(i)=par(i)/a
      end do
      nsp=ns/2
      print 900,nsp,nl,(ntsl(i),i=1,nl)
  900 format(' SETKA:   sum : ntsl(nl) =',i5,' lines = ',i4/
     *       10i5/10i5)
      do i=1,ids
        dolm=ddolgs*(i-1)
        do j=1,its
            tetm=dtets*(j-1)
          call ggmraw(1,dolg,tet,dolm,tetm)
          gkoor(1,j,i)=tet
          gkoor(2,j,i)=dolg
        end do
      end do
      return
      end
