      subroutine lambet(ne,i1,i2,det,ht,tt,co2,cn2,co,ch,che,
     *           tn,vnq,cim,cio,cih,cihe,vio,vih,vihe,ti,te,col,
     *           vdu,vdv,beta,lam,vio1,vih1,vihe1,al,ga,
     *           cio1,cih1,cihe1,dolm)
C     LAM - diffusion coefficients for atomic ions
C           (NE=1 - O+, NE=2 - H+, NE=3 - HE+)
C     or conductivities ions (NE=4) and electrons (NE=5)
C     BETA - transport coefficients atomic ions or
C     convection heatings:
C       GA-AL*FU-KAP*D/DS(-LAM*D/DS(FU)-BETA*FU)=0
C                        or
C       GA-AL*FU-KAP*D/DS(-LAM*D/DS(FU))-BETA*D/DS(FU)=0
      USE mo_bas_gsm, ONLY: nv0,re
      real lam,n2,nu4,nu5,nu,la

      dimension ht(nv0),tt(nv0),co2(nv0),cn2(nv0),co(nv0),ch(nv0),che(nv0),
     *          tn(nv0),vnq(nv0),cim(nv0),cio(nv0),cih(nv0),col(nv0),
     *          cihe(nv0),vio(nv0),vih(nv0),vihe(nv0),ti(nv0),te(nv0),
     *          vdu(nv0),vdv(nv0),beta(nv0),lam(nv0),vio1(nv0),vih1(nv0),
     *          vihe1(nv0),al(nv0),ga(nv0),cio1(nv0),cih1(nv0),cihe1(nv0)

C    DET - time step 
      eht=1./det
      do 7 m=i1,i2
        alt=ht(m)
        tet=tt(m)
        tei=ti(m)
        ten=tn(m)
        o2 = 0.
        n2 = 0.
        o  = 0.
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
c  Transport coefficients
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
C                    branch O+ (NE=1)
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
c                      for HE+ (NE=3)
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
C  Calculation of the flux of thermal electrons for the convective transfer term 
C  in the heat balance equation
        if(ne.gt.3)w=oi*vio(m)+hi*vih(m)+hei*vihe(m)
c       if(ne.gt.3)w=0.
C     Convecting heating coefficient BETA for 
c     heat balance equation of ions
        if(ne.eq.4)beta(m)=w/ces
C                           ... & electrons
        if(ne.eq.5)beta(m)=w/ce
c  Classic electron condactivity and diffusion
        lam(m)=la(ne,o2,n2,o,h,he,oi,hi,hei,ce,bm,rps,tei,tee,alt)
c  Classic electron condactivity end

    7 continue

      return
      end
      
C           Calculation of the ion-neutral friction term 
C           in the equations of motion of atomic ions     
      function rik(ne,ti,tn,co2,cn2,co,ch,che,alt)
      real nu
      data c11/6.67e-10/,c12/6.87e-10/,c14/1.29e-10/,c15/1.32e-10/,
     *c21/3.21e-9/,c22/3.39e-9/,c23/2.37e-9/,c25/1.06e-9/,
     *c31/2.11e-9/,c32/2.19e-9/,c33/1.09e-9/,c34/4.74e-10/,
     *c24/3.35e-10/,c35/8.77e-11/,c13/5.59e-11/
      goto(1,2,3),ne
    1 continue
C       O+ - neutrals (NE=1)
        a=0.
        if(alt.le.1.e8)a=c11*co2+c12*cn2
        b=0.
        if(alt.le.5.e8)b=nu(1,ti,tn)*co*c13
        rik=a+b+c14*ch+c15*che
        return
    2 continue
C       H+ - neutrals (NE=2)
        a=0.
        if(alt.le.1.e8)a=c21*co2+c22*cn2
        b=0.
        if(alt.le.5.e8)b=c23*co
        rik=a+b+nu(2,ti,tn)*ch*c24+c25*che
        return
    3 continue
C       HE+ - neutrals (NE=3)
        a=0.
        if(alt.le.1.e8)a=c31*co2+c32*cn2
        b=0.
        if(alt.le.5.e8)b=c33*co
        rik=a+b+c34*ch+nu(1,ti,tn)*che*c35
        return
      end

C             Calculation of the ion-ion friction term 
C             in the equations of motion of atomic ions
      function sik(ne,cio,cih,cihe,nu5,nu4)
      real nu4,nu5
      data c1/2.47e-6/,c2/3.95e-5/,c3/1.256/,c4/1.9e-1/,
     *c5/7.6e-1/,c6/3.14e-1/
      goto(1,2,3,4,5,6),ne
    1 continue
C        O+ - H+ (NE=1)
        sik=c1*cih
        return
    2 continue
C        O+ - HE+ (NE=2)
        sik=c4*nu5*cihe
        return
    3 continue
C        H+ - O+ (NE=3)
        sik=c2*cio
        return
    4 continue
C        H+ - HE+ (NE=4)
        sik=c3*nu4*cihe
        return
    5 continue
C        HE+ - O+ (NE=5)
        sik=c5*nu5*cio
        return
    6 continue
C        HE+ - H+ (NE=6)
        sik=c6*nu4*cih
        return
      end

