      subroutine betnv(bet,bet0,tv,tef)
c . . . calculation vibrational T N2
      dimension a(8),b(8)
      data a/3.39e-15,2.33e-14,3.02e-14,-2.74e-14,
     *      -3.84e-15,1.6e-14,-2.3e-14,  2.77e-14/
     *    ,b/3.72e-13,3.09e-11,1.92e-10,2.9e-10,
     *       5.85e-11,1.59e-10,1.19e-10,1.36e-10/
     *    ,e1/3353./
      tev=e1/tv
      x0=1.-exp(-tev)
      ski=0.
      do 1 i=1,8
       x=x0*exp(-tev*i)
c    . . . quasitrinor distribution
c      x=x0*exp(-i*tev+i*(i-1)*20.6/tef)
       ak=a(i)*tef+b(i)
       ski=ski+x*ak
    1 continue
      bet=bet0*x0+ski
      return
      end
