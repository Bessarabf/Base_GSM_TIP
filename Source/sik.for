      function sik(ne,cio,cih,cihe,nu5,nu4)
      real nu4,nu5
      data c1/2.47e-6/,c2/3.95e-5/,c3/1.256/,c4/1.9e-1/,
     *c5/7.6e-1/,c6/3.14e-1/
      goto(1,2,3,4,5,6),ne
    1 continue
        sik=c1*cih
        return
    2 continue
        sik=c4*nu5*cihe
        return
    3 continue
        sik=c2*cio
        return
    4 continue
        sik=c3*nu4*cihe
        return
    5 continue
        sik=c5*nu5*cio
        return
    6 continue
        sik=c6*nu4*cih
        return
      end

