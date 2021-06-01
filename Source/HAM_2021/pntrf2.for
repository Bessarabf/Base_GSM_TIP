      function pntrf2(alt,co,tn,te)
      dimension eps(3),a(3),b(3),c(3),e(3),dx(3),ex(3)
      data c1/5.34e-10/,c2/2.71e-10/,c3/2.63e-8/,c4/7.e-5/
      data eps/0.02,0.028,0.08/,a/7.883e-6,9.466e-6,1.037e-8/
      data b/1.021,0.8458,1.633/,c/1.009,0.9444,1.466/
      data e/228.,326.,98./
c Pavlov, Berrington
c      d=5+exp(-326.6/tn)+3.*exp(-227.7/tn)
c      f=8.132e-9/d
c      pntrf2=c3/3.4e-12*f*co/(te*tn)
c Pavlov, Berrington
c Dalgarno, Degges
c      pntrf2=c3*co*(1.-c4*te)/tn
c Dalgarno, Degges
c a la Bailey
c       pntrf2=0.5*c3*co*(1.-c4*te)/tn
c a la Bailey
c Bailey, Balan
       t1=tn
       t0=tn
       z=5.+3.*exp(-228/t1)+exp(-326/t0)
       dx(1)=exp(-228./t1)
       dx(2)=exp(-326./t0)
       dx(3)=exp(-326./t0)
       ex(1)=exp(-228./te)
       ex(2)=exp(-326./te)
       ex(3)=exp(-98./te-228./t1)
       s=0.
       do i=1,3
         f=(1.+b(i))*dx(i)+(e(i)/te+1.+b(i))*ex(i)
         f=5.91e-9*(te-tn)*f+eps(i)*(ex(i)-dx(i))
         s=s+a(i)*c(i)*te**(b(i)-0.5)*f
       end do
       pntrf2=c3/3.4e-12*8.629e-6*co/z*s
c Bailey, Balan
      return
      end
