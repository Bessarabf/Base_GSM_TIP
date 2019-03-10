      subroutine conducn(bb,cm,co2,cn2,co,tn,sp,sh)
      data e/1.60219e-20/,oe/1.76e7/,oi/3.09e2/,
     *ci1/4.23e-10/,ci2/4.28e-10/,ci3/2.58e-10/,
     *ce1/1.82e-10/,ce11/3.6e-2/,ce2/2.33e-11/,ce21/1.21e-4/,
     *ce3/2.8e-10/
      ee=e/bb*cm
      ome=oe*bb
      omi=oi*bb
      b=sqrt(tn)
      fi=ci1*co2+ci2*cn2+ci3*co
      fe=ce1*(1.+ce11* b)*b*co2+ce2*(1.-ce21*tn)*tn*cn2+ce3*b*co
      bu=fi/omi
      cu=fe/ome
      a=1./(1.+bu*bu)
      b=1./(1.+cu*cu)
      sp=ee*(bu*a+cu*b)
      sh=ee*(a-b)
      return
      end
