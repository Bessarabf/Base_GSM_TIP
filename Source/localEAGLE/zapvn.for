c
      subroutine zapvn(i,kpart,par,vdr,nr,ks,vn,ins)
      dimension par(nr),vdr(ks),vn(ins)
      i1=i*kpart+1
      i2=i*2+1
      vn(1)=par(i1)
      vn(2)=par(i1+3)
      vn(3)=vdr(i2)
      vn(4)=vdr(i2+1)
      vn(5)=par(i1+7)
      vn(6)=par(i1+6)
      return
      end

