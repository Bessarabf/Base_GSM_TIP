      function piqj(alt,co2,cn2,co,ch,che,tn,vnq,vnu,vnv,
     *cio,cih,cihe,vio,vih,vihe,ti,vdu,vdv)
      real nu
      data c1/2.44e-8/,c2/3.17e-13/,c3/1.01e-8/
      f=nu(3,ti,tn)
      vs=(vdu-vnu)**2+(vdv-vnv)**2
c     piqj1=vs*cio*pqji(1,alt,co2,cn2,co,ch,che,tn,ti)
c     piqj2=vs*cih*pqji(2,alt,co2,cn2,co,ch,che,tn,ti)
c     piqj3=vs*cihe*pqji(3,alt,co2,cn2,co,ch,che,tn,ti)
c     piqj5=0.
c     piqj4=0.
c     piqj6=0.
      piqj1=((vio-vnq)**2+vs)*cio*pqji(1,alt,co2,cn2,co,ch,che,tn,ti)
      piqj2=((vih-vnq)**2+vs)*cih*pqji(2,alt,co2,cn2,co,ch,che,tn,ti)
      piqj3=((vihe-vnq)**2+vs)*cihe*pqji(3,alt,co2,cn2,co,ch,che,tn,ti)
      piqj5=cio*cihe*(vio-vihe)**2*c1*f
      piqj4=cio*cih*(vio-vih)**2*c2
      piqj6=cih*cihe*(vih-vihe)**2*c3*f
      piqj=piqj1+piqj2+piqj3+piqj4+piqj5+piqj6
      return
      end
