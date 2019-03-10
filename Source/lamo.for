      real function lamo(ne,tn,vnq,vnu,vnv,ti,vio,
     *vdu,vdv,te)
c    *vdu,vdv)
      data c11/2.82e-11/,c12/7.74e-12/,c13/1.073e-12/,c14/5.17e-14/,
     *c15/9.65e-16/,c21/1.533e-12/,c22/5.92e-13/,c23/8.6e-14/,
     *c24/2.73e-12/,c25/1.155e-12/,c26/1.483e-13/,c31/4.3e-11/,
     *c1/6.67e-1/,c2/4.28e-8/,c3/6.36e-1/,c4/4.08e-8/
      goto(1,1,4),ne
    1 continue
        vs=(vio-vnq)**2+(vdu-vnu)**2+(vdv-vnv)**2
        goto(2,3),ne
    2   continue
          t=c2*vs+c1*(ti-tn)+tn
          a=t/300.
          as=a*a
          aq=as*a
          af=aq*a
          lamo=c11-c12*a+c13*as-c14*aq+c15*af
          return
    3   continue
          t=c4*vs+c3*(ti-tn)+tn
c
          tv=te
c
c         tv=tn
c         tv=tn*1.15
cb        tv=tn*1.4
c         tv=tn*1.5
c         tv=tn*1.25
c         tv=tn*1.75
          a=t/300.
          as=a*a
c         if(t.le.1700.)lamo=c21-c22*a+c23*as
c         if(t.gt.1700.)lamo=c24-c25*a+c26*as
          if(t.le.1700.)bet0=c21-c22*a+c23*as
          if(t.gt.1700.)bet0=c24-c25*a+c26*as
          call betnv(bet,bet0,tv,t)
          lamo=bet
          return
    4 continue
        a=sqrt(tn+ti*6.25e-2)
        lamo=c31*a
        return
      end