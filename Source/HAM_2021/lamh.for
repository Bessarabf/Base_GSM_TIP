      real function lamh(tn,ti)
      data c1/3.8e-11/
      a=sqrt(ti+tn*6.25e-2)
      lamh=c1*a
      return
      end
