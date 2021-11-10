      subroutine zatm(shi,dol,ut,fzatm)
      data pi/3.1415926/
      data re/6370./, x0/3300./
      data utna/34200./, utko/45360./
      fshi(xx) = 363.011 - 28.529*xx - 5.248e03/(xx**3) -
     -   7.605*xx*xx
      fdol(xx) = 9.78e03 - 3.672e03*xx - 6.513e04/(xx**3) +
     + 371.882*xx*xx
**
** shi,dol - географические  широта и долгота узла в градусах
** текущее время  UT в секундах
** utna,utko-время начала 9h30m и конца 12h36m затмения в сек
** на расстоянии х0 от центра затмения нет
**
      doli = dol
      if (ut.lt.utna)  then
         fzatm = 1.
         return
      end if
      if (ut.gt.utko)  then
         fzatm = 1.
         return
      end if
**
      xxx = ut/10000.
      fi0 = fdol(xxx)
      sh0 = fshi(xxx)
**
      if (fi0.gt.180.) fi0=fi0-360.
      if (doli.gt.180.) doli=doli-360.
          dfi = (doli-fi0)*pi/180.
          dshi = (shi-sh0)*pi/180.
      ro = cos(shi*pi/180.)*dfi
**
***   x1 = re*sqrt(ro*ro + dshi*dshi)
      x2 = re*acos(cos(ro)*cos(dshi))
**
**    print *,x1,x2
                x = x2
      if (x.ge.x0) then
      fzatm = 1.
      return
      else if(x.lt.x0)  then
      y = -0.03*x + 100.
      fzatm = 1.-y/100.
      return
      end if
      end
