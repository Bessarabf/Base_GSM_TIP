      subroutine difte(dd,alt,co2,cn2,co,ch,che,tes,tn,ti,
     *           cim,cio,cih,cihe,dgte)
      dte=dd*tes
      te1=tes-dte
      te2=tes+dte
            p1=pnte(alt,co2,cn2,co,ch,che,te1)
            p2=pntrf(alt,co2,cn2,co,tn,te1)
            p3=pmte(cim,te1)
            p4=petd12(co2,cn2,tn,te1,tn)
            p5=petd3(co,tn,te1)
            p6=pite(cio,cih,cihe,te1)
      gte1=-(p1+p2+p6)*(te1-tn)-p3*(te1-ti)+p4+p5
            p1=pnte(alt,co2,cn2,co,ch,che,te2)
            p2=pntrf(alt,co2,cn2,co,tn,te2)
            p3=pmte(cim,te2)
            p4=petd12(co2,cn2,tn,te2,tn)
            p5=petd3(co,tn,te2)
            p6=pite(cio,cih,cihe,te2)
      gte2=-(p1+p2+p6)*(te2-tn)-p3*(te2-ti)+p4+p5
      dgte=(gte2-gte1)/(2.*dte)
      return
      end


