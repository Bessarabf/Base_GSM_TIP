      subroutine njuan(qq,uu,oo,tt,h16)
      double precision x(4),q,u,c,s,r,p,v,a,b,apb,amb,o,t,rmv
      data re/6371.02e5/
      q=qq
      u=uu
      c=dble(1.)/dble(3.)
      s=u*u*.125/q
      r=s*s*.25
      p=dble(1.)/dble(27.)+r
      r=p+r
      v=dabs(s)*dsqrt(p)
      a=(r+v)**c
      rmv=r-v
      b=dabs(rmv)
      b=rmv/b*b**c
      v=a+b
      apb=v*.5
      amb=(a-b)*.5
      p=amb*dsqrt(dble(3.))
      r=c-apb
      a=datan2(p,r)*.5
      p=dsqrt(v+c)
      b=(amb*amb*3.+r*r)**.25*2.*dcos(a)
      apb=-p
      amb=-b
      x(1)=p+b
      x(2)=p+amb
      x(3)=apb+b
      x(4)=apb+amb
      h0=h16*1.e-5
      d0=0.d1
      do1i=1,4
        if(dabs(x(i)).gt.1.d0)goto1
          if(x(i).le.d0.and.s.gt.d0.or.x(i).ge.d0.and.s.lt.d0)goto1
            a=x(i)
            goto2
    1 continue
    2 continue
      s=dble(1.)-a*a
      p=dsqrt(s)
      t=datan2(p,a)
      o=u/s
      tt=t
      oo=re*(1./o-1.)
      return
      end

