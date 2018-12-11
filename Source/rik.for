      function rik(ne,ti,tn,co2,cn2,co,ch,che,alt)
      real nu
      data c11/6.67e-10/,c12/6.87e-10/,c14/1.29e-10/,c15/1.32e-10/,
     *c21/3.21e-9/,c22/3.39e-9/,c23/2.37e-9/,c25/1.06e-9/,
     *c31/2.11e-9/,c32/2.19e-9/,c33/1.09e-9/,c34/4.74e-10/,
     *c24/3.35e-10/,c35/8.77e-11/,c13/5.59e-11/
      goto(1,2,3),ne
    1 continue
        a=0.
        if(alt.le.1.e8)a=c11*co2+c12*cn2
        b=0.
        if(alt.le.5.e8)b=nu(1,ti,tn)*co*c13
        rik=a+b+c14*ch+c15*che
        return
    2 continue
        a=0.
        if(alt.le.1.e8)a=c21*co2+c22*cn2
        b=0.
        if(alt.le.5.e8)b=c23*co
        rik=a+b+nu(2,ti,tn)*ch*c24+c25*che
        return
    3 continue
        a=0.
        if(alt.le.1.e8)a=c31*co2+c32*cn2
        b=0.
        if(alt.le.5.e8)b=c33*co
        rik=a+b+c34*ch+nu(1,ti,tn)*che*c35
        return
      end
