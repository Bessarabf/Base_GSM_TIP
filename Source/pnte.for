      function pnte(alt,co2,cn2,co,ch,che,te)
      data c1/1.92e-14/,c2/9.13e-16/,c3/6.24e-15/,c4/4.94e-12/,
     *c5/1.26e-13/,c6/1.21e-4/,c7/3.6e-2/,c8/1.35e-4/
      ts=sqrt(te)
      pnte2=0.
      pnte1=0.
      pnte3=0.
      if(alt.gt.1.e8)goto1
        pnte2=c2*cn2*(1.-c6*te)*te
        pnte1=c3*co2*(1.+c7*ts)*ts
    1 continue
      if(alt.le.5.e8)pnte3=c1*co*ts
      pnte4=c4*ch*(1.-c8*te)*ts
      pnte5=c5*che*ts
      pnte=pnte1+pnte2+pnte3+pnte4+pnte5
      return
      end
