      subroutine cxy(ce,hsl,nol,nu,nvn,tsl,pga,pun)
      dimension ce(*),hsl(nvn,*),nu(*),tsl(nvn,*)
      dimension x(140),y(140),z(140)
      data re/6371.2/,pd/1.74532925e-2/,sq3/1.7320508/
      ns1=19
      nsn=nu(nol)-18
      r=re+hsl(ns1,nol)
      a=tsl(ns1,nol)*pd
      s=sin(a)
      scx=r/(s*s*2.)
      c=cos(a)
      s=sqrt(c*c*3.+1.)
        d1=s
        rk1=r*r*r
      q=rk1/s
      scy=q*.75
      i=0
      do j=ns1,nsn
        i=i+1
        tet=tsl(j,nol)*pd
        c=cos(tet)
        ca=sqrt(c*c*3.+1.)
        r=re+hsl(j,nol)
        a=c*sq3+ca
        x(i)=(alog(a)/sq3+ca*c)*scx
          rk=r*r*r
        q=rk/ca
        a=1.-scy/q
        y(i)=ce(j)/sqrt(a)
        z(i)=ce(j)*d1/ca*rk/rk1*1.e5
      end do
      pun=unta(1,i,x,y)*5.2e-7
      pga=unta(1,i,x,z)*1.04e-11
      end
