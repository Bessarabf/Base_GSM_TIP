!     Non elastic heat lost electrons on neutrals
      function pntrf(alt,co2,cn2,co,tn,te)
      data c1/5.34e-10/,c2/2.71e-10/,c3/2.63e-8/,c4/7.e-5/
!     lost on rotate exitation Ž2
      pntrf1=0.
!     lost on rotate exitation N2
      pntrf2=0.
      if(alt.gt.1.e8)goto1
        ts=1./sqrt(te)
        pntrf1=c1*co2*ts
        pntrf2=c2*cn2*ts
    1 continue
!     lost on thin structure exitate Ž
      pntrf3=c3*co*(1.-c4*te)/tn
      pntrf=pntrf1+pntrf2+pntrf3
      return
      end
