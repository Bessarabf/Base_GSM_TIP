      subroutine backfv(ne,i1,i2,ht,tt,cim,cio,cih,cihe,ti,te,
     *                  bet,gam,ga,al,lam)
      real lam
      dimension ht(*),tt(*),cim(*),cio(*),cih(*),ti(*),
     *          cihe(*),bet(*),gam(*),ga(*),al(*),lam(*),te(*)
      data pi/3.1415926/
      i3=i2-1
c     print*,' subrotine backfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
      do4mn=i1,i3
        m=i3-mn+i1
        mp=m+1
        mm=m-1
        if(m.eq.i1)mm=m
        alt=ht(m)
        bm=bdip(alt,tt(m))
        if(ne.gt.3)ce=cio(m)+cih(m)+cihe(m)
        if(ne.eq.5.and.alt.le.1.e8)ce=ce+cim(m)
        if(ne.gt.3)bm=bm/ce
        dst=ds(ht(mm),alt,ht(mp),tt(mm),tt(mp))
        dst=dst/(bm+bm)
        fl=lam(mp)-dst*ga(m)
c       if(ne.gt.3.and.mn.eq.i3)fl=0.
        fu=gam(m)-bet(m)*fl
c       if(m.ne.i1)lam(m)=fl+dst*al(m)*fu
        if(m.ne.i1)lam(m)=fl*(1.-bet(m)*dst*al(m))+dst*al(m)*gam(m)
        goto(1,2,3,5,6),ne
    1   continue
          cio(m)=fu
c         if(cio(m).lt.1.e-3)cio(m)=1.e-3
          goto4
    2   continue
          cih(m)=fu
c         if(cih(m).lt.1.e-3)cih(m)=1.e-3
          goto4
    3   continue
          cihe(m)=fu
c         if(cihe(m).lt.1.e-3)cihe(m)=1.e-3
          goto4
    5   continue
          ti(m)=fu
          goto4
    6   continue
          te(m)=fu
    4 continue
c     print*,' subrotine backfv'
c     print*,' ne=',ne
c     print*,' i1,i2=',i1,i2
c     do i=i1,i2
c       write(*,10)ht(i)*1.e-5,tt(i)*180./pi,cio(i),cih(i),lam(i)
c  10   format(f8.0,f8.2,3(1pe12.2))
c     end do
c     print*,'backfv output'
      return
      end

