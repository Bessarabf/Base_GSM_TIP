!     Lost electron heat on exitation stage 0
!             
      function petd3(co,tn,te)
      data c1/1.21e-8/,c2/3000./,c3/-22713./,c4/2.4e4/,
     *     c5/3.e-1/,c6/1500./,c7/1.947e-5/,c8/4000./

      a=(te-tn)/(te*tn)
      if(a.lt.0) a=0
      b=exp(c3*a)-1.
      c=(te-c2)/(c2*te)
      d=te-c6
      e=te-c8
      f=c4+c5*d-c7*d*e
      a=exp(f*c)
      petd3=c1*co*a*b
      return
      end
