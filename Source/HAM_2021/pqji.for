      function pqji(ne,alt,co2,cn2,co,ch,che,tn,ti)
      real nu
      data c11/5.7e-17/,c12/5.61e-17/,c13/3.58e-18/,c14/9.72e-19/,
     *c15/3.39e-18/,c21/2.50e-17/,c22/2.62e-17/,c23/1.79e-17/,
     *c24/1.34e-18/,c25/6.83e-18/,c31/6.01e-17/,c32/6.14e-17/,
     *c33/2.79e-17/,c34/3.04e-18/,c35/1.41e-18/
      goto(1,2,3),ne
    1 continue
        a=0.
        b=0.
        if(alt.le.1.e8)a=c11*co2+c12*cn2
        if(alt.le.5.e8)b=c13*nu(1,ti,tn)*co
        pqji=a+b+c14*ch+c15*che
        return
    2 continue
        a=0.
        b=0.
        if(alt.le.1.e8)a=c21*co2+c22*cn2
        if(alt.le.5.e8)b=c23*co
        pqji=a+b+c24*nu(2,ti,tn)*ch+c25*che
        return
    3 continue
        a=0.
        b=0.
        if(alt.le.1.e8)a=c31*co2+c32*cn2
        if(alt.le.5.e8)b=c33*co
        pqji=a+b+c34*ch+c35*nu(1,ti,tn)*che
        return
      end

