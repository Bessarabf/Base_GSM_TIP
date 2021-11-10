      function pntrf1(alt,co2,cn2,tn,te)
      data c1/5.34e-10/,c2/2.71e-10/,c3/2.63e-8/,c4/7.e-5/
      pntrf2=0.
      pntrf3=0.
      if(alt.gt.1.e8)goto1
        ts=1./sqrt(te)
        pntrf2=c1*co2*ts
        pntrf3=c2*cn2*ts
    1 continue
      pntrf1=pntrf2+pntrf3
      return
      end
