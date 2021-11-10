!    electrons heat losses for vibrational excitation of neutrals
!
      function petd12(co2,cn2,tn,te,tv)
      data c11/5.76e-9/,c12/3.902e3/,c13/4.38e2/,c14/4.56e-4/,
     *c15/2400./,c16/700./,c17/-3000./,c21/2.31e-8/,c22/2000./,
     *c23/1.06e4/,c24/7.51e3/,c25/1.1e-3/,c26/1800./,c27/3300./,
     *c28/1.233/,c29/1000./,c30/2.056e-4/,c31/4000./
c     tv=te
cb    tv=tn*1.4
c     tv=tn*1.5
c     tv=tn*1.25
c     tv=tn
c     tv=tn*1.15
c     tv=tn*1.75
      a=(te-tn)/(te*tn)
      b=exp(c17*a)-1.
      c=c13*tanh(c14*(te-c15))
      d=exp((c12+c)*(te-c16)/(c16*te))
c   Vibrational excitation of O2 
      ptd1=c11*co2*d*b
      b=te-c29
      c=te-c31
      d=te-c26
      f=c23+c24*tanh(c25*d)
      g=c27+c28*b-c30*b*c
      a=(te-tv)/(te*tv)
      b=exp(-g*a)-1.
      c=exp(f*(te-c22)/(c22*te))
c . . . and N2
      ptd2=c21*cn2*c*b
      petd12=ptd1+ptd2
      return
      end
