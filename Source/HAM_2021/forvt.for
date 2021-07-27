      subroutine forvt(ne,i1,i2,ht,tt,cim,cio,cih,cihe,beta,bet,
     *                 gamm,ga,al,lam,betap,gap,
     *                 alp,lamp,betam,gam,alm,lamm,temp)
      real lam,lamp,lamm

      dimension ht(*),tt(*),cim(*),cio(*),cih(*),cihe(*),
     *          beta(*),bet(*),gamm(*),ga(*),al(*),lam(*),
     *          betap(*),betam(*),gap(*),gam(*),alp(*),alm(*),
     *          lamp(*),lamm(*),temp(*)
c     data dte/.5/
      data dte/.1/
      i3=i1+1
      i4=i2-1
      do1m=i3,i4
c       dte=temp(m)*5.e-5
        mp=m+1
        mm=m-1
        dtem=dte*2.*temp(m)
        dtemp=dte*2.*temp(mp)
        dtemm=dte*2.*temp(mm)
        alt=ht(m)
        altp=ht(mp)
        altm=ht(mm)
        tet=tt(m)
        tetp=tt(mp)
        tetm=tt(mm)
        bg=bdip(alt,tet)
        ce=cio(m)+cih(m)+cihe(m)
        if(ne.eq.5.and.alt.le.1.e8)ce=ce+cim(m)
        dst=ds(altm,alt,altp,tetm,tetp)
        bg=bg/(ce*dst)
        hp=(altp+alt)*.5
        dst=ds(alt,hp,altp,tet,tetp)
        gp=dst
        bm=beta(m)
        bp=abs(bm)
        a=lam(m)
        d=-(bm-bp)*.5/dst
        f=d*(temp(mp)-temp(m))
        e=bg*(a+lam(mp))/dst
        b=e+d
        f=f+e*(temp(mp)-temp(m))
        hm=(altm+alt)*.5
        dst=ds(altm,hm,alt,tetm,tet)
        gm=dst
        d=(bm+bp)*.5/dst
        f=f-d*(temp(m)-temp(mm))
        e=bg*(a+lam(mm))/dst
        f=f-e*(temp(m)-temp(mm))
        a=e+d
        c=a+b+al(m)
        c=c+(alp(m)-alm(m))/dtem*temp(m)-(gap(m)-gam(m))/dtem
        c=c+bg*(lamp(m)-lamm(m))*(temp(m)-temp(mm))/(gm*dtem)
        c=c-bg*(lamp(m)-lamm(m))*(temp(mp)-temp(m))/(gp*dtem)
        a=a-bg*(lamp(mm)-lamm(mm))/dtemm*(temp(m)-temp(mm))/gm
        b=b+bg*(lamp(mp)-lamm(mp))/dtemp*(temp(mp)-temp(m))/gp
        hp=betap(m)
        hm=betam(m)
        d=abs(hp)
        e=abs(hm)
        c=c+(hp-d-hm+e)*.5/dtem*(temp(mp)-temp(m))/gp
        c=c+(hp+d-hm-e)*.5/dtem*(temp(m)-temp(mm))/gm
        f=f+ga(m)-al(m)*temp(m)
        c=1./(c-a*bet(mm))
        bet(m)=b*c
        gamm(m)=(f+a*gamm(mm))*c
    1 continue
      return
      end


