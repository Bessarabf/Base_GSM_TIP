      function pitn(alt,co2,cn2,co,ch,che,tn,cio,
     *cih,cihe,ti)
      real nu
      data c11/4.44e-10/,c12/5.00e-10/,c13/5.6e-11/,c14/2.42e-10/,
     *c15/2.11e-10/,c21/1.95e-10/,c22/2.34e-10/,c23/2.79e-10/,
     *c25/4.26e-10/,c31/4.68e-10/,c32/5.47e-10/,c33/4.35e-10/,
     *c34/7.58e-10/,c24/3.35e-10/,c35/8.75e-11/
      a=0.
      b=0.
      c=0.
      if(alt.gt.1.e8)goto1
        a=c11*co2+c12*cn2
        b=c21*co2+c22*cn2
        c=c31*co2+c32*cn2
    1 continue
      d=0.
      e=0.
      f=0.
      if(alt.gt.5.e8)goto2
        d=c13*nu(1,ti,tn)*co
        e=c23*co
        f=c33*co
    2 continue
c
      pitn1=a+d+c14*ch+c15*che
c
      pitn2=b+e+c24*nu(2,ti,tn)*ch+c25*che
c
      pitn3=c+f+c34*ch+c35*nu(1,ti,tn)*che
      pitn=cio*pitn1+cih*pitn2+cihe*pitn3
      return
      end

