c     ver 18.04.14 add nv to interface
      subroutine alga(ne,nx,i1,i2,dt,ht,tt,co2,cn2,co,ch,che,tn,
     *                vnq,vnu,vnv,cim,cio,cih,cihe,vio,vih,vihe,
     *                ti,te,vdu,vdv,cio1,cih1,cihe1,ti1,te1,qo,qsm,
     *		      al,ga,qomt,qmax,iqo,mass,NV)

      dimension ht(NV),tt(NV),co2(NV),cn2(NV),co(NV),ch(NV),
     *          che(NV),tn(NV),vnq(NV),vnu(NV),vnv(NV),cim(NV),
     *          cio(NV),cih(NV),cihe(NV),vio(NV),vih(NV),vihe(NV),
     *          ti(NV),te(NV),vdu(NV),vdv(NV),
     *          cio1(NV),cih1(NV),cihe1(NV),ti1(NV),te1(NV),
     *          qo(NV),al(NV),ga(NV),qsm(NV),
     *          ,mass(30)
      dimension hs(23),hn(23),yos(23),yon(23),tns(23),tnn(23),tes(23),
     *          tenor(23),cnes(23),cnen(23),alfa(23),tvn(23),tvs(23)
      real la1,la2,la3,la4,n2,lamh,lamo,mi,la5,la6
      data la5/1.e-9/,la6/1.4e-9/,e0/10./,qt0/2.0e14/
	!


      if(nx.ge.45) then
        do m=1,23        !!!!!!!!!!!!!!!!!!!!!!!!!!
          m1=nx-m+1
          hs(m)=ht(m)*1.e-5
!
          hn(m)=ht(m1)*1.e-5
          yos(m)=co(m)
          yon(m)=co(m1)
          tns(m)=tn(m)
          tnn(m)=tn(m1)
          tes(m)=te(m)
          tenor(m)=te(m1)
          cnes(m)=cim(m)+cio(m)+cih(m)+cihe(m)
          cnen(m)=cim(m1)+cio(m1)+cih(m1)+cihe(m1)
        end do
        call altv(hs,yos,tns,tes,cnes,alfa,tvs,23)
        call altv(hn,yon,tnn,tenor,cnen,alfa,tvn,23)
c       print*,' nx=',nx,' i1=',i1
c       print*,' Tvs'
c       print*,tvs
c       print*,' Tvn'
c       print*,tvn
      end if
      if(ne.ne.5) goto 17
        qt=0.
        if(nx.eq.NV) goto 17
          i3=(nx+1)/2
          if(i3.le.18) goto 17
            call ntnb(nx,ht,tt,cio,cih,cihe,cim,b1,snt,snb)
            edd=alog10(snb)-11.
            ed=1.2
            if(edd.gt.0.)ed=1.3*edd**2.2+ed
c           ed=sqrt(4./3.*2.6e-12*snb)
            a=ed/e0
            qt=qt0*(1.-(1.+a)*exp(-a))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           qt=qt*1.5
ccc         qt=qt*2.0
c           qt=qt*3.0
c           if(nx.ge.63)then
c           if(nx.ge.79)then
c           if(nx.ge.93)then
c           if(nx.ge.109)then
c             qt=qt*(7.*nx-745.)/18.
c             qt=qt*(nx-107.)/2.
c линейный рост в 3 раза от 60 град широты до 70 град широты
c             qt=qt*(nx-100.)/9.
c
c линейный рост в 5 раз от 60 град широты до 70 град широты
c             qt=qt*(2.*nx-209.)/9.
c
c линейный рост в 3,5 раза от 60 град широты до 70 град широты
c             qt=qt*(5.*nx-509.)/36.
c
c             qt=qt*(9.*nx-476.)/64.
c             qt=qt*(9.*nx-663.)/48.
c             qt=qt*(9.*nx-803.)/34.
c             qt=qt*(9.*nx-476.)/128.
c           end if
c
c поток фотоэлектронов на 70 град широты увеличен в 6 раз
c           if(nx.eq.127)qt=qt*4.
c!!         if(nx.eq.127)qt=qt*6.
c!!         if(nx.eq.127)qt=qt*8.
c
c
c поток фотоэлектронов с 65 град широты увеличен в 4 раза
c           if(nx.ge.117)qt=qt*4.
c
c
c поток фотоэлектронов с 65 град широты увеличен в 4.5 раза
c           if(nx.ge.117)qt=qt*4.5
c
c
c поток фотоэлектронов с 65 град широты увеличен в 6 раз
c           if(nx.ge.117)qt=qt*6.0
c
c           qt=qt*(19.*nx+99.)/118.
c           qt=qt0*(1.-(1.+a)*exp(-a))*10.
c           a=(qo(1)+qo(nx))*.5
c           if(a.lt.10.)qt=qt*.1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
            qt=qt/(snt*b1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           qt=qt/qomt
            qt=qt/qmax
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
   17 continue
      eht=1./dt
      do16m=i1,i2
        alt=ht(m)
        h=ch(m)
        vud=vdu(m)
        tet=tt(m)
        dv=dvn(alt,tet,vud)
        if(ne.eq.3)goto1
          tei=ti(m)
          ten=tn(m)
          voi=vio(m)
          if(alt.le.5.e8)o=co(m)
    1   continue
        if(ne.eq.1)goto2
          he=che(m)
    2   continue
        if(ne.eq.2)goto3
          if(alt.gt.1.e8)goto3
            n2=cn2(m)
            o2=co2(m)
    3   continue
        if(ne.eq.3.or.ne.eq.5)goto4
          vqn=vnq(m)
          vun=vnu(m)
          vvn=vnv(m)
          vvd=vdv(m)
    4   continue
        if(ne.eq.1.or.ne.eq.3)goto5
          oi=cio(m)
    5   continue
        if(ne.eq.2.or.ne.eq.3)goto6
          hi=cih(m)
    6   continue
        if(ne.eq.1.and.alt.le.1.e8.or.ne.eq.5.and.alt.le.1.e8)oq=qo(m)
        if(ne.eq.5.and.alt.le.1.e8)qmol=qsm(m)
        goto(7,7,8,13,13),ne
    7   continue
          if(alt.le.5.e8)la4=lamh(ten,tei)*o
c
          tv=ten
          if(nx.ge.45)then
            if(m.ge.1.and.m.le.23)tv=tvs(m)
            if(m.ge.nx-22.and.m.le.nx)tv=tvn(nx-m+1)
          end if
          la3=lamo(3,ten,vqn,vun,vvn,tei,voi,vud,vvd,tv)
c         la3=lamo(3,ten,vqn,vun,vvn,tei,voi,vud,vvd,te(m))
c         la3=lamo(3,ten,vqn,vun,vvn,tei,voi,vud,vvd)
c
          goto9
    8   continue
c         ga=n(he+)*eht
c         al=la5*n(o2)+la6*n(n2)+dv+eht,
          ga(m)=cihe1(m)*eht
c
c         ga(m)=0.
c
          if(dv.lt.0.)ga(m)=ga(m)-cihe1(m)*dv
          a=eht
          if(dv.ge.0.)a=dv+a
          if(alt.le.1.e8)a=a+la5*o2+la6*n2
          al(m)=a
          goto16
    9   continue
        goto(10,12),ne
   10   continue
c         ga=oq+la4*n(o)*n(h+)+n(o+)*eht
c         al=la1*n(o2)+la2*n(n2)+la3*n(h)+dv+eht,
          g=cio1(m)*eht
c
c         g=0.
c
          if(dv.lt.0.)g=g-cio1(m)*dv
          a=la3*h+eht
          if(dv.ge.0.)a=a+dv
          if(alt.gt.1.e8)goto11
c
          tv=ten
          if(nx.ge.45)then
            if(m.ge.1.and.m.le.23)tv=tvs(m)
            if(m.ge.nx-22.and.m.le.nx)tv=tvn(nx-m+1)
          end if
            la1=lamo(1,ten,vqn,vun,vvn,tei,voi,vud,vvd,tv)
            la2=lamo(2,ten,vqn,vun,vvn,tei,voi,vud,vvd,tv)
c           la1=lamo(1,ten,vqn,vun,vvn,tei,voi,vud,vvd,te(m))
c           la2=lamo(2,ten,vqn,vun,vvn,tei,voi,vud,vvd,te(m))
c           la1=lamo(1,ten,vqn,vun,vvn,tei,voi,vud,vvd)
c           la2=lamo(2,ten,vqn,vun,vvn,tei,voi,vud,vvd)
c
            g=g+oq+la5*o2*cihe(m)
            a=a+la1*o2+la2*n2
   11     continue
          if(alt.le.5.e8)g=g+la4*hi
          ga(m)=g
          al(m)=a
          goto16
   12   continue
c         ga=la3*n(h)*n(o+)+n(h+)*eht
c         al=la4*n(o)+dv+eht
          ga(m)=la3*h*oi+cih1(m)*eht
c
c         ga(m)=la3*h*oi
c
          if(dv.lt.0.)ga(m)=ga(m)-cih1(m)*dv
          a=eht
          if(dv.ge.0.)a=a+dv
          if(alt.le.5.e8)a=a+la4
          al(m)=a
          goto16
   13   continue
          mm=m-1
          mp=m+1
          if(m.eq.i1)mm=m
          if(m.eq.i2)mp=m
          hei=cihe(m)
          if(alt.le.1.e8)mi=cim(m)
          tee=te(m)
	!!!!!!!
	    if(tee.lt.ten) tee=ten	!!!!!
          pte=pite(oi,hi,hei,tee)
          ces=oi+hi+hei
          ce=ces
          if(alt.le.1.e8)ce=ce+mi
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
c     ga=(ptn*tn+ce*pte*te+pqj)/ces+ti*eht
c     al=(ptn+ce*pte+ra)/ces+dvt+eht
c     dvt=2/3*dv,
c     ces=n(o+)+n(h+)+n(he+),
c     ce=ces+n(m+)
            ptn=pitn(alt,o2,n2,o,h,he,ten,oi,hi,hei,tei)

            pqj=piqj(alt,o2,n2,o,h,he,ten,vqn,vun,vvn,
     *      oi,hi,hei,voi,vhi,vhei,tei,vud,vvd)
c           pqj=0.
            ces=1./ces
            b=ce*pte
            rao=oi*(viop*bmp-viom*bmm)*dst
            rah=hi*(vihp*bmp-vihm*bmm)*dst
            rahe=hei*(vihep*bmp-vihem*bmm)*dst
            ra=6.67e-1*bm*(rao+rah+rahe)
c           ra=0.
            g=(ptn*ten+pqj+b*tee)*ces+ti1(m)*eht
            if(ra.lt.0.)g=g-ti1(m)*ra*ces
            if(dvt.lt.0.)g=g-ti1(m)*dvt
c           a=(ptn+ra+b)*ces+eht
            a=(ptn+b)*ces+eht
            if(ra.ge.0.)a=a+ra*ces
            if(dvt.ge.0.)a=a+dvt
            ga(m)=g
            al(m)=a
            goto16
   15   continue
c     ga=(ptn*tn+pte*ti+ptd-ptd1-ptd2)+te*eht
c     al=ptn+pte+ra+dvt+eht
          q=0.
          if(alt.le.1.e8)q=oq
          ptn=pnte(alt,o2,n2,o,h,he,tee)
c старый вариант
c          if(alt.le.5.e8)ptn=ptn+pntrf(alt,o2,n2,o,ten,tee)
c старый вариант
c новый вариант
          if(alt.le.5.e8)ptn=ptn+pntrf1(alt,o2,n2,ten,tee)
c кроме Bailey, Balan
c          if(alt.le.5.e8)ptn=ptn+pntrf2(alt,o,ten,tee)
c кроме Bailey, Balan
c новый вариант
c
          tv=ten
          if(nx.ge.45)then
            if(m.ge.1.and.m.le.23)tv=tvs(m)
            if(m.ge.nx-22.and.m.le.nx)tv=tvn(nx-m+1)
          end if
c
          if(alt.le.1.e8)ptd1=petd12(o2,n2,ten,tee,tv)
c
c         if(alt.le.1.e8)ptd1=petd12(o2,n2,ten,tee)
c
          if(alt.le.5.e8)ptd2=petd3(o,ten,tee)
          ciom=cio(mm)
          cihm=cih(mm)
          cihem=cihe(mm)
          ciop=cio(mp)
          cihp=cih(mp)
          cihep=cihe(mp)
          cem=ciom+cihm+cihem
          cep=ciop+cihp+cihep
          if(hm.le.1.e8)cem=cem+cim(mm)
          if(hp.le.1.e8)cep=cep+cim(mp)
          vem=(ciom*viom+cihm*vihm+cihem*vihem)*bmm/cem
          vep=(ciop*viop+cihp*vihp+cihep*vihep)*bmp/cep
          ra=6.67e-1*bm*(vep-vem)*dst
c         ra=0.
          if(alt.le.1.e8)q=q+qmol
c          ptd=pgfel(q,ce,alt)
c
c дополнительный источник локального нагрева с 65 град широты
c         if(nx.ge.117)ptd=ptd*2.
c
c
c дополнительный источник локального нагрева с 70 град широты
ccc       if(nx.eq.127)ptd=ptd*2.3
c         if(nx.gt.127)ptd=ptd*3.
c
          ptd=pgfkr(q,ce,alt,o2,n2,o,h,he)
c          ptd=pgfkhaz(q,ce,alt,o2,n2,o,h,he)
c         if(alt.gt.7.e7)ptd=ptd+qt*bm*(qo(1)+qo(nx))*5.e-3/8.
          if(alt.gt.7.e7)ptd=ptd+qt*bm*(qo(iqo)+qo(nx-iqo+1))*.5
          g=ptn*ten+pte*tei+ptd+te1(m)*eht
c новый вариант
c Bailey, Balan
          g=g-pntrf2(alt,o,ten,tee)
c Bailey, Balan
c новый вариант
          if(ra.lt.0.)g=g-te1(m)*ra
          if(dvt.lt.0.)g=g-te1(m)*dvt
          if(alt.le.1.e8)g=g+ptd1
          if(alt.le.5.e8)g=g+ptd2
          ga(m)=g
c         a=ptn+pte+eht+ra
          a=ptn+pte+eht
          if(ra.ge.0.)a=a+ra
          if(dvt.ge.0.)a=a+dvt
          al(m)=a
   16 continue
     
      return
      end
