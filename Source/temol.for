      subroutine temol(par,pari,rads,mass,kpars,nh,its,dts,ntr,ins)
      dimension par(kpars,nh,its)
     *          ,mass(30),rads(nh),pari(ins,nh,its)
      real la1,la2,lamo
      data alf/1.26e-4/,na/1/,dd/0.10/
      ch=0.
      che=0.
      cih=0.
      cihe=0.
      ntr2=ntr+2
!      print *,'nh=',nh,'ntr2=',ntr2
      do 1 ig=2,its-1
        do 2 i=1,ntr2
          alt=rads(i)
          tn=par(7,i,ig)
          co2=par(1,i,ig)
          cn2=par(2,i,ig)
          co=par(3,i,ig)
          qo=par(16,i,ig)
          cim=par(6,i,ig)
          te0=par(9,i,ig)
          ti=par(8,i,ig)
          qs=par(13,i,ig)+par(14,i,ig)+par(15,i,ig)
          qs=qs+qo
          pqs=qs/cim
          pqc=alf*cim
c
          tia=tn
          if(i.ge.16)tia=pari(6,i,ig)
          la1=lamo(1,tn,0.,0.,0.,tia,0.,0.,0.,tn)
          la2=lamo(2,tn,0.,0.,0.,tia,0.,0.,0.,tn)
c         la1=lamo(1,tn,0.,0.,0.,tia,0.,0.,0.,te0)
c         la2=lamo(2,tn,0.,0.,0.,tia,0.,0.,0.,te0)
c         la1=lamo(1,tn,0.,0.,0.,tn,0.,0.,0.)
c         la2=lamo(2,tn,0.,0.,0.,tn,0.,0.,0.)
c
          cio=qo/(la1*co2+la2*cn2)
          ce=cim+cio
          qe=pgfkr(qs,ce,alt,co2,cn2,co,ch,che)
          tes=te0
        if(tes.lt.0.) print22,tes,i,ig
  22    format(' tes=',e10.3,' i= ',i4,' ig=',i4)
          m=1
    4     if(m.gt.mass(16))go to 3
            p1=pnte(alt,co2,cn2,co,ch,che,tes)
            p2=pntrf(alt,co2,cn2,co,tn,tes)
            p3=pmte(cim,tes)
            p4=petd12(co2,cn2,tn,tes,tn)
            p5=petd3(co,tn,tes)
            p6=pite(cio,cih,cihe,tes)
            f1=(p1+p2+p6)*tn+p3*ti+p4+p5+pqc
            f2=p1+p2+p3+p6+pqs
            call difte(dd,alt,co2,cn2,co,ch,che,tes,tn,ti,
     *                    cim,cio,cih,cihe,dgte)
c           tes=(te0+(f1+qe)*dts)/(1.+f2*dts)
            tes=(te0+(f1+qe-tes*dgte)*dts)/(1.+(f2-dgte)*dts)
            if(tes.le.tn) tes=tn
            m=m+1
            go to 4
    3     continue
          par(9,i,ig)=tes
    2   continue
        ntr3=ntr+3
        do 5 i=ntr3,nh
          par(9,i,ig)=pari(5,i,ig)
    5   continue
    1 continue
      return
      end


