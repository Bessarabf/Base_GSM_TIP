      subroutine lambet(ne,i1,i2,det,ht,tt,co2,cn2,co,ch,che,
     *           tn,vnq,cim,cio,cih,cihe,vio,vih,vihe,ti,te,col,
     *           vdu,vdv,beta,lam,vio1,vih1,vihe1,al,ga,
     *           cio1,cih1,cihe1,dolm)

      real lam,n2,nu4,nu5,nu,la

      dimension ht(*),tt(*),co2(*),cn2(*),co(*),ch(*),che(*),
     *          tn(*),vnq(*),cim(*),cio(*),cih(*),col(*),
     *          cihe(*),vio(*),vih(*),vihe(*),ti(*),te(*),
     *          vdu(*),vdv(*),beta(*),lam(*),vio1(*),vih1(*),
     *          vihe1(*),al(*),ga(*),cio1(*),cih1(*),cihe1(*)
      data re/6371.02e5/
c     print*,'lambet inpout'
c     print*,'ne=',ne,'i1=',i1,'i2=',i2,'dolm=',dolm
      eht=1./det
      do 7 m=i1,i2
        alt=ht(m)
        tet=tt(m)
        tei=ti(m)
        ten=tn(m)
        if(alt.le.1.e8) o2=co2(m)
        if(alt.le.1.e8) n2=cn2(m)
        if(alt.le.5.e8) o=co(m)
        h=ch(m)
        he=che(m)
        oi=cio(m)
        hi=cih(m)
        hei=cihe(m)
        tee=te(m)
        bm=bdip(alt,tet)
        ce=oi+hi+hei
        ces=ce
        if(alt.le.1.e8)ce=ce+cim(m)
        goto(1,1,1,6,6),ne
    1   continue
c  Classic electron condactivity
          mm=m-1
          mp=m+1
          if(m.eq.i1)mm=m
          if(m.eq.i2)mp=m
          hp=ht(mp)
          hm=ht(mm)
          dst=ds(hm,alt,hp,tt(mm),tt(mp))
          dst=1./dst
c  Classic electron condactivity end
          vud=vdu(m)
          gs=gst(alt,tet,vdv(m),vud,dolm)
          ci=cos(tet)
          ci=sin(tet)/sqrt(1.+3.*ci*ci)
          r=re+alt
          gs=gs-(2.*vud*oml(dolm)+omr(tet,dolm)*ci*
     *    (2.*vdv(m)-r*omt(tet,dolm)))
          dti=(ti(mp)-ti(mm))*dst
c  Classic electron condactivity
          dte=(te(mp)-te(mm))*dst
c  Classic electron condactivity end
          dnm=0.
          if(alt.lt.1.e8.and.hp.le.1.e8.and.hm.le.1.e8)
     *    dnm=(cim(mp)-cim(mm))*dst
          if(alt.le.1.e8.and.hp.gt.1.e8)
     *    dnm=(cim(m)-cim(mm))*dst*2.
          if(alt.le.1.e8.and.hm.gt.1.e8)
     *    dnm=(cim(mp)-cim(m))*dst*2.
          dnv=dvn(alt,tet,vud)
cccc
          rin=0.
cccc      rin=eht
cccc
c         rin=dnv
c         rin=dnv+eht
cccc
cccc      if(dnv.ge.0.)rin=rin+dnv*.5
cccc
          goto(2,3,4),ne
    2     continue
            r=rik(1,tei,ten,o2,n2,o,h,he,alt)
!	print*,'lambet ',ne,i1,i2,tei,ten
            nu5=nu(3,tei,ten)
            s1=sik(1,oi,hi,hei,nu5,nu4)
            s2=sik(2,oi,hi,hei,nu5,nu4)
            dnk=(cih(mp)-cih(mm))*dst
            dnl=(cihe(mp)-cihe(mm))*dst
            fob=1./oi
c           if(vio(m).ge.0.)dv=(vio(m)-vio(mm))*dst*2.
c           if(vio(m).lt.0.)dv=(vio(mp)-vio(m))*dst*2.
c           rin=rin+dv
c           rin=rin+(ga(m)*fob-al(m))+(1.-cio1(m)*fob)*eht
cccc
cccc        gs=gs+vio1(m)*eht
cccc        if(dnv.lt.0.)gs=gs-.5*dnv*vio(m)
cccc
            goto5
    3     continue
c
            r=rik(2,tei,ten,o2,n2,o,h,he,alt)
            nu4=nu(3,tei,ten)
            s1=sik(3,oi,hi,hei,nu5,nu4)
            s2=sik(4,oi,hi,hei,nu5,nu4)
            dnk=(cio(mp)-cio(mm))*dst
            dnl=(cihe(mp)-cihe(mm))*dst
            fob=1./hi
c           if(vih(m).ge.0.)dv=(vih(m)-vih(mm))*dst*2.
c           if(vih(m).lt.0.)dv=(vih(mp)-vih(m))*dst*2.
c           rin=rin+dv
c           rin=rin+(ga(m)*fob-al(m))+(1.-cih1(m)*fob)*eht
cccc
cccc        gs=gs+vih1(m)*eht
cccc        if(dnv.lt.0.)gs=gs-.5*dnv*vih(m)
cccc
            goto5
    4     continue
c
            r=rik(3,tei,ten,o2,n2,o,h,he,alt)
            nu5=nu(3,tei,ten)
            nu4=nu5
            s1=sik(5,oi,hi,hei,nu5,nu4)
            s2=sik(6,oi,hi,hei,nu5,nu4)
            dnk=(cio(mp)-cio(mm))*dst
            dnl=(cih(mp)-cih(mm))*dst
            fob=1./hei
c           if(vihe(m).ge.0.)dv=(vihe(m)-vihe(mm))*dst*2.
c           if(vihe(m).lt.0.)dv=(vihe(mp)-vihe(m))*dst*2.
c           rin=rin+dv
c           rin=rin+(ga(m)*fob-al(m))+(1.-cihe1(m)*fob)*eht
cccc
cccc        gs=gs+vihe1(m)*eht
cccc        if(dnv.lt.0.)gs=gs-.5*dnv*vihe(m)
cccc
    5     continue
          rps=1./((r+s1+s2+rin)*bm)
c         rps=1./bm
c         col(m) =r+s1+s2+rin
          r=r*vnq(m)
          dn=dnk+dnl
          if(alt.le.1.e8)dn=dn+dnm
          dt=tee/ce*dn+dti+dte
          r=(r+gs)*rps
          dt=dt*rps
          beta(m)=be(ne,dt,r,rps,vio(m),vih(m),vihe(m),s1,s2)
    6   continue
        if(ne.gt.3)w=oi*vio(m)+hi*vih(m)+hei*vihe(m)
c       if(ne.gt.3)w=0.
        if(ne.eq.4)beta(m)=w/ces
        if(ne.eq.5)beta(m)=w/ce
c  Classic electron condactivity
        lam(m)=la(ne,o2,n2,o,h,he,oi,hi,hei,ce,bm,rps,tei,tee,alt)
c  Classic electron condactivity end

    7 continue
c     print*,'lambet output'
      return
      end
