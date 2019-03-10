      real function la(ne,co2,cn2,co,ch,che,cio,cih,cihe,
     *ce,bm,rps,ti,te,alt)
      real la1,la2,la3,la4,la5
      data c11/5.2e6/,c21/8.31e7/,c31/2.08e7/,c41/1.84e8/,
     *c42/9.18e7/,c43/3.72e8/,c51/5.94e9/,c52/7.08e-12/,c53/2.54e-13/,
     *c54/9.02e-13/,c55/1.09e-16/,c56/1.09e-11/,c57/1.76e-10/,
     *c58/2.4e-14/,c59/1.8e-11/,s/2.5/
      goto(1,2,3,4,5),ne
    1 continue
c
        t=ti+cio/ce*te
        la=c11*t*rps
        return
    2 continue
c
        t=ti+cih/ce*te
        la=c21*t*rps
        return
    3 continue
c
        t=ti+cihe/ce*te
        la=c31*t*rps
        return
    4 continue
c
        t=ti**s/bm
        c=(c41*cihe+c42*cio+c43*cih)/ce
        la=t*c
        return
    5 continue
c
        t=te**s*c51/bm
        la1=0.
        la2=0.
        la3=0.
        if(alt.gt.1.e8)goto6
          sq=sqrt(te)
          la1=(c52+c53*sq)*co2
          la2=(c54-c55*te)*sq*cn2
    6   continue
        if(alt.le.5.e8)la3=c56*co
        la4=(c57-c58*te)*ch
        la5=c59*che
        c=1.+te*te/ce*(la1+la2+la3+la4+la5)
c       c=1.+sqrt(2.)+te*te/ce*(la1+la2+la3+la4+la5)
        la=t/c
        return
      end
