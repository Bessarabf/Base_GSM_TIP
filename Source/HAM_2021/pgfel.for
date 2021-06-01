      function pgfel(q,ce,alt)
      data c1/7.74e3/,c2/1.e-5/,c3/120./,c4/1.875e-2/,
     *c5/7.74e4/
      if(alt.le.6.e7)eps=c1*(1.+(alt*c2-c3)*c4)
      if(alt.gt.6.e7)eps=c5
      pgfel=0.
      if(alt.le.1.e8) pgfel=eps*q/ce
      return
      end
