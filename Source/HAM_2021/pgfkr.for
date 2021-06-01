      function pgfkr(q,ce,alt,o2,n2,o,h,he)
      real n2
      data c1/1.238e3/
      pgfkr=0.
      if(alt.gt.1.e8)goto1
        cn=o2+n2+o+h+he
        cn=ce/cn
        a=alog10(cn)
        if(a.gt.-3.)a=-3.
        if(a.lt.-8.)a=-8.
        pgfkr=c1/ce*(9.+a)**2*q
    1 continue
      return
      end
