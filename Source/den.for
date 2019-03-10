      function den(denl,difn,tl1,tz1,alpha,weight,eddy,eddy1,
     *             omega,ratio,zh,z1,h1,r2,z2,h2)
      a=2.8e1/(2.895e1-weight)
      a1=1.e0/a
      aexp=exp(a)
      b=alog(denl*difn*(tl1/tz1)**(alpha+1.e0))
      c=alog(eddy*eddy1*(tl1/tz1))
      ab=a*b
      ab=abs(ab)
      ac=a*c
      ac=abs(ac)
      if  (ab.gt.1.70e2.or.ac.gt.1.70e2)  goto  1
      d=0.e0
      den=exp(d)*(aexp**(b-d-10.)+aexp**(c-d-10.))**a1*
     *    aexp**(10.*a1)
      goto  2
  1   d=amax1(ab,ac)
      d=d/a
      d=aint(d)
      d=abs(d)
      den=exp(d)*(aexp**(b-d)+aexp**(c-d))**a1
  2   continue
      den=den*exp(alog(omega*ratio)/(1.e0+exp((zh-z1)/h1)))
      if  (zh.le.2.40e2)  den=den*exp(r2/(1.e0+exp((zh-z2)/h2)))
      if (den.lt.0.01) then
        print *,' p/p den. den=0.'
        stop
      end if
      return
      end

                                                                        
