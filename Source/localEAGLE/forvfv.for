      subroutine forvfv(ne,i1,i2,ht,tt,cim,cio,cih,cihe,
     *alf,beta,bet,gam,ga,al,lam,vio,vih,vihe,col)
      real lam,lamp
      dimension ht(*),tt(*),cim(*),cio(*),cih(*),
     *cihe(*),bet(*),beta(*),gam(*),ga(*),al(*),
     *lam(*),vio(*),vih(*),vihe(*),col(*)
      data pi/3.14159265/,sigm/.93/
c     data pi/3.14159265/,sigm/.87/
c     print*,' subrotine forvfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
      i3=i1+1
      altm=ht(i1)
      altp=ht(i3)
      tetm=tt(i1)
      tetp=tt(i3)
      dstm=ds(altm,altm,altp,tetm,tetp)
      bm=bdip(altm,tetm)
      if(ne.gt.3)ce=cio(i1)+cih(i1)+cihe(i1)
      if(ne.eq.5.and.altm.le.1.e8)ce=ce+cim(i1)
      if(ne.gt.3)bm=bm/ce
      dstm=dstm/(bm+bm)
      i4=i2-1
      do1m=i1,i4
        g=ga(m)
        mp=m+1
        mm=m-1
        mn=mp+1
        if(m.eq.i1)mm=m
        if(m.eq.i4)mn=mp
        lamp=(lam(m)+lam(mp))*.5
cc      coll=(col(m)+col(mp))*.5
        alt=(altm+altp)*.5
        dst=1./ds(altm,alt,altp,tetm,tetp)
        alt=lamp*dst
        betp=(beta(m)+beta(mp))*.5
        ab=abs(betp)
        bp=(betp+ab)*.5
        lamp=alt+bp
        bm=(betp-ab)*.5
        betp=alt-bm
        dstp=ds(altm,altp,ht(mn),tetm,tt(mn))
        bgm=bdip(altm,tetm)
        bgmp=bdip(altp,tetp)
        bgmm=bdip(ht(mm),tt(mm))
        if(ne.gt.3)ce=cio(mp)+cih(mp)+cihe(mp)
        if(ne.eq.5.and.altp.le.1.e8)ce=ce+cim(mp)
        if(ne.gt.3)bgmp=bgmp/ce
        dstp=dstp/(bgmp+bgmp)
        a=al(mp)
        bmp=dstp
        bmm=dstm
        b=bet(m)
c       write(*,11)ht(m)*1.e-5,tt(m)*180./pi,lamp
c       write(*,11)ht(m)*1.e-5,tt(m)*180./pi,alt
c  11   format(f8.0,f8.2,1pe15.5)
        if(lamp.lt.1.)then
          alt=1.+b*betp
cc        alt=coll+b*betp
          dst=1./(lamp+bmp*alt*a)
          alf=dst*lamp
          bet(mp)=dst*alt
c         bmm=betp*(gam(m)+bmm*b*g/sigm)
c         goto(2,3,4),ne
c   2     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp*b*(bmm*al(m)*cio(m)+
c    *        .5*(cio(mp)*vio(mp)/bgmp-cio(mm)*vio(mm)/bgmm))+
c    *        .5*(cio(m)*vio(m)/bgm+cio(mp)*vio(mp)/bgmp)+
c    *        lamp*cio(mp)-betp*cio(m))
c           goto5
c   3     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp*b*(bmm*al(m)*cih(m)+
c    *        .5*(cih(mp)*vih(mp)/bgmp-cih(mm)*vih(mm)/bgmm))+
c    *        .5*(cih(m)*vih(m)/bgm+cih(mp)*vih(mp)/bgmp)+
c    *        lamp*cih(mp)-betp*cih(m))
c           goto5
c   4     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp*b*(bmm*al(m)*cihe(m)+
c    *        .5*(cihe(mp)*vihe(mp)/bgmp-cihe(mm)*vihe(mm)/bgmm))+
c    *        .5*(cihe(m)*vihe(m)/bgm+cihe(mp)*vihe(mp)/bgmp)+
c    *        lamp*cihe(mp)-betp*cihe(m))
c   5     continue
cc        gam(mp)=dst*bmm
          gam(mp)=dst*betp*(gam(m)+bmm*b*g)
        else
          alt=1./lamp+betp/lamp*b
cc        alt=coll/lamp+betp/lamp*b
          dst=1./(1.+bmp*alt*a)
          alf=dst
          bet(mp)=dst*alt
cc        bmm=betp/lamp*(gam(m)+bmm*b*g/sigm)
c         goto(6,7,8),ne
c   6     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp/lamp*b*(bmm*al(m)*cio(m)+
c    *        .5*(cio(mp)*vio(mp)/bgmp-cio(mm)*vio(mm)/bgmm))+
c    *        .5*(cio(m)*vio(m)/bgm+cio(mp)*vio(mp)/bgmp)/lamp+
c    *        cio(mp)-betp/lamp*cio(m))
c           goto9
c   7     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp/lamp*b*(bmm*al(m)*cih(m)+
c    *        .5*(cih(mp)*vih(mp)/bgmp-cih(mm)*vih(mm)/bgmm))+
c    *        .5*(cih(m)*vih(m)/bgm+cih(mp)*vih(mp)/bgmp)/lamp+
c    *        cih(mp)-betp/lamp*cih(m))
c           goto9
c   8     continue
c           bmm=bmm-(1.-sigm)/sigm*(betp/lamp*b*(bmm*al(m)*cihe(m)+
c    *        .5*(cihe(mp)*vihe(mp)/bgmp-cihe(mm)*vihe(mm)/bgmm))+
c    *        .5*(cihe(m)*vihe(m)/bgm+cihe(mp)*vihe(mp)/bgmp)/lamp+
c    *        cihe(mp)-betp/lamp*cihe(m))
c   9     continue
cc        gam(mp)=dst*bmm
          gam(mp)=dst*betp/lamp*(gam(m)+bmm*b*g)
        end if
        altm=altp
        altp=ht(mn)
        tetm=tetp
        tetp=tt(mn)
        dstm=dstp
    1 continue
c     print*,' subrotine forvfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
c     print*,' alf=',alf
c     do i=i1,i2
c       write(*,10)ht(i)*1.e-5,tt(i)*180./pi,bet(i),gam(i)
c  10   format(f8.0,f8.2,2(1pe15.5))
c     end do
      return
      end
