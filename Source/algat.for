      subroutine algat(ne,nx,m,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                    vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *                    ti,te,vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *                    al,ga,qomt,qmax,iqo,mass,NV)

      dimension ht(NV),tt(NV),co2(NV),cn2(NV),co(NV),ch(NV),
     *          che(NV),tn(NV),vnq(NV),vnu(NV),vnv(NV),cim(NV),
     *          cio(NV),cih(NV),cihe(NV),vio(NV),vih(NV),vihe(NV),
     *          vdu(NV),vdv(NV),cio1(NV),cih1(NV),cihe1(NV),
     *          ti1(NV),te1(NV),qo(NV),qsm(NV),
     *          mass(30)
!	* ,al(NV),ga(NV)
      real la1,la2,la3,la4,n2,lamh,lamo,mi,la5,la6
      data la5/1.e-9/,la6/1.4e-9/,e0/10./

      eht=1./dt
      alt=ht(m)
      h=ch(m)
      vud=vdu(m)
      tet=tt(m)
      dv=dvn(alt,tet,vud)
      tei=ti
      ten=tn(m)
      voi=vio(m)
!	print*,'algat m=',m,nv,nx,i1,i2,co
!	pause

      o=co(m)
      he=che(m)
      n2=cn2(m)
      o2=co2(m)
      if(ne.eq.5)goto4
        vqn=vnq(m)
        vun=vnu(m)
        vvn=vnv(m)
        vvd=vdv(m)
    4 continue
      oi=cio(m)
      hi=cih(m)
      if(ne.eq.5)oq=qo(m)
      if(ne.eq.5)qmol=qsm(m)
      mm=m-1
      mp=m+1
      if(m.eq.i1)mm=m
      if(m.eq.i2)mp=m
      hei=cihe(m)
      mi=cim(m)
      tee=te
      if (tee.lt.ten) TEE=TEN !!!!!
      pte=pite(oi,hi,hei,tee)
      ces=oi+hi+hei
      ce=ces
      ce=ce+mi
      bm=bdip(alt,tet)
      hp=ht(mp)
      hm=ht(mm)
      tp=tt(mp)
      tm=tt(mm)
      vhi=vih(m)
      vhei=vihe(m)
      dst=1./ds(hm,alt,hp,tm,tp)
      bmp=1./bdip(hp,tp)
      bmm=1./bdip(hm,tm)
      dvt=6.67e-1*dv
      w=oi*voi+hi*vhi+hei*vhei
      viop=vio(mp)
      viom=vio(mm)
      vihp=vih(mp)
      vihm=vih(mm)
      vihep=vihe(mp)
      vihem=vihe(mm)
      if(ne.eq.5)goto15
        ptn=pitn(alt,o2,n2,o,h,he,ten,oi,hi,hei,tei)
        pqj=piqj(alt,o2,n2,o,h,he,ten,vqn,vun,vvn,
     *  oi,hi,hei,voi,vhi,vhei,tei,vud,vvd)
        ces=1./ces
        b=ce*pte
        rao=oi*(viop*bmp-viom*bmm)*dst
        rah=hi*(vihp*bmp-vihm*bmm)*dst
        rahe=hei*(vihep*bmp-vihem*bmm)*dst
        ra=6.67e-1*bm*(rao+rah+rahe)
        g=(ptn*ten+pqj+b*tee)*ces+ti1(m)*eht
        if(ra.lt.0.)g=g-ti1(m)*ra*ces
        if(dvt.lt.0.)g=g-ti1(m)*dvt
        a=(ptn+b)*ces+eht
        if(ra.ge.0.)a=a+ra*ces
        if(dvt.ge.0.)a=a+dvt
        ga=g
        al=a
        goto16
   15 continue
        q=oq
        ptn=pnte(alt,o2,n2,o,h,he,tee)
c старый вариант
c        ptn=ptn+pntrf(alt,o2,n2,o,ten,tee)
c старый вариант
c новый вариант
        ptn=ptn+pntrf1(alt,o2,n2,ten,tee)
c кроме Bailey, Balan
c        ptn=ptn+pntrf2(alt,o,ten,tee)
c кроме Bailey, Balan
c новый вариант
        ptd1=petd12(o2,n2,ten,tee,ten)
        ptd2=petd3(o,ten,tee)
        ciom=cio(mm)
        cihm=cih(mm)
        cihem=cihe(mm)
        ciop=cio(mp)
        cihp=cih(mp)
        cihep=cihe(mp)
        cem=ciom+cihm+cihem
        cep=ciop+cihp+cihep
        cem=cem+cim(mm)
        cep=cep+cim(mp)
        vem=(ciom*viom+cihm*vihm+cihem*vihem)*bmm/cem
        vep=(ciop*viop+cihp*vihp+cihep*vihep)*bmp/cep
        ra=6.67e-1*bm*(vep-vem)*dst
        q=q+qmol
c        ptd=pgfel(q,ce,alt)
c
c дополнительный источник локального нагрева с 65 град широты
c       if(nx.ge.117)ptd=ptd*2.
c
c
c дополнительный источник локального нагрева с 70 град широты
ccc     if(nx.eq.127)ptd=ptd*2.3
c       if(nx.gt.127)ptd=ptd*3.
c
c        ptd=pgfkhaz(q,ce,alt,o2,n2,o,h,he)
        ptd=pgfkr(q,ce,alt,o2,n2,o,h,he)
        g=ptn*ten+pte*tei+ptd+te1(m)*eht
c новый вариант
c Bailey, Balan
        g=g-pntrf2(alt,o,ten,tee)
c Bailey, Balan
c новый вариант
        if(ra.lt.0.)g=g-te1(m)*ra
        if(dvt.lt.0.)g=g-te1(m)*dvt
        g=g+ptd1
        g=g+ptd2
        ga=g
        a=ptn+pte+eht
        if(ra.ge.0.)a=a+ra
        if(dvt.ge.0.)a=a+dvt
        al=a
   16 continue
      return
      end
