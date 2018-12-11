      subroutine int2pef(x1,x2,xi,vn1,vn2,plm)
      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
      plm=w1*vn2+w2*vn1
      return
      end
                                           
