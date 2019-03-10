      function erf(x)
      double precision y,t,u,p,a1,a2,a3,a4,a5,v,w
      data p,a1,a2,a3,a4,a5/.3275911,.254829592,-.284496736,
     *1.421413741,-1.453152027,1.061405429/
      y=x
      t=1.d0/(1.d0+p*y)
      u=dexp(-y*y)
      v=t*t
      w=a1*t+a2*v
      v=v*t
      w=w+a3*v
      v=v*t
      w=w+a4*v
      v=v*t
      w=w+a5*v
      erf=1.d0-w*u
      return
      end
