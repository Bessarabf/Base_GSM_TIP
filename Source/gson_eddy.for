      subroutine gson_eddy(an1,an31,an2,an3,vr,an6,q,
     *                  rp,r,g,lv,dt,eddyco,ro,n1,n2,n,vi,vj)
c . . .
      dimension an1(n1,n2,n),an2(n1,n2,n),
     *          an3(n1,n2,n),an6(n1,n2,n),an31(n1,n2,n),
     *          vr(n1,n2,n),q(n1,n2,n),rp(n),r(n),g(n),eddyco(n,n1),
     *          ro(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n),
     *          dp(3,3)
      allocatable aal(:),bbc(:),dd31(:),dd32(:),ql(:),
     *          fp(:),ap(:),bp(:),ep(:),ed(:),b2(:),
     *          b3(:),cp(:),dpp(:),pfu(:),ppo(:)
      data am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
      data bk/1.38e-16/,pi/3.14159/,re/6.371e08/
      
	
      allocate (aal(n),bbc(n),dd31(n),dd32(n),ql(n),
     *          fp(n),ap(n),bp(n),ep(n),ed(n),b2(n),
     *          b3(n),cp(n),dpp(n),pfu(n),ppo(n))


      if(n.lt.lv) print *,'mistake on gsozkg. n can''t low lv'
      np=lv-1
      n11=n1-1
      n21=n2
      dtet=pi/(n1-1)
      dfi=2.*pi/n2
      gam=1.e-20
c*******
      do 40 i=2,n11
       tetl=dtet*(i-1)
       do 401 j=1,n21
        jp=j+1
        jl=j-1
        if(j.eq.1)jl=n2
        if(j.eq.n21)jp=1
c
        do 1 k=1,lv
            ams=(32.*an1(i,j,k)+an2(i,j,k)*28.+16.*an3(i,j,k))
     *         /(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))*1.67e-24
            cono2=an1(i,j,k)
            conn2=an2(i,j,k)
            cono=an3(i,j,k)
            tempp=an6(i,j,k)
            rog=ro(i,j,k)
c
c           call kmold(dp,cono2,conn2,cono,tempp)
            call kmoldn(dp,cono2,conn2,cono,tempp,rog)
            hh=bk*tempp/g(k)
            hcp=hh/ams
            ho=hh/am3
            hoo=hh/am1
            hnn=hh/am2
            kv=k+1
            kn=k-1
            if(k.eq.lv)kv=lv
            if(k.eq.1)kn=1
            gtem=(an6(i,j,kv)-an6(i,j,kn))/(r(kv)-r(kn))
            aal(k)=eddyco(k,i)+dp(3,3)
            bbc(k)=eddyco(k,i)*(1./hcp+1./tempp*gtem)+
     *      dp(3,3)*(1./ho+1./tempp*gtem)
c           bbc(k)=bbc(k)-vr(i,j,k)*0.
cc          b2(k)=dp(3,1)*(1./hoo+1./tempp*gtem)
cc            b3(k)=dp(3,2)*(1./hnn+1./tempp*gtem)
cc            b2(k)=b2(k)/an3(i,j,k)*an1(i,j,k)
cc            b3(k)=b3(k)/an3(i,j,k)*an2(i,j,k)
cc            goo=(an1(i,j,kv)-an1(i,j,kn))/(r(kv)-r(kn))
cc   *        *    dp(3,1)/an3(i,j,k)
cc            gnn=(an2(i,j,kv)-an2(i,j,kn))/(r(kv)-r(kn))
cc   *        *    dp(3,2)/an3(i,j,k)
cc            b2(k)=b2(k)+goo
cc            b3(k)=b3(k)+gnn
              rt=an6(i,j,k)*(r(kv)-r(kn))
              ht=(an1(i,j,kv)*an6(i,j,kv)-
     *            an1(i,j,kn)*an6(i,j,kn))/rt
              hb=an1(i,j,k)/hoo
                hs=ht+hb
c               hm=hs/hb
                hm=hs/ht
cc            if(abs(hm).le.0.1)hs=0.
                b2(k)=dp(3,1)*hs/an3(i,j,k)

              ht=(an2(i,j,kv)*an6(i,j,kv)-
     *            an2(i,j,kn)*an6(i,j,kn))/rt
              hb=an2(i,j,k)/hnn
                hs=ht+hb
c               hm=hs/hb
                hm=hs/ht
cc            if(abs(hm).le.0.1)hs=0.
                b3(k)=dp(3,2)*hs/an3(i,j,k)
              dd31(k)=dp(3,1)
              dd32(k)=dp(3,2)
 1        continue
c*******
          do 3 k=2,np
            alf=9.9e-34*exp(470./an6(i,j,k))
            bet=1.1e-34*exp(510./an6(i,j,k))
            pp=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
            if(q(i,j,k).ne.0.)bet=0.
            qq=2.*q(i,j,k)*an1(i,j,k)
            ql(k)=2.*alf*an3(i,j,k)*pp+bet*an1(i,j,k)*pp
     *      +2.*gam*an3(i,j,k)
c*******
            aa=0.5/dfi/(r(k)+re)/sin(tetl)
            dd=aa*(vj(i,jp,k)-vj(i,jl,k))
            aa=0.5/dtet/(r(k)+re)/sin(tetl)
c            dd=dd+aa*(vi(i+1,j,k)*sin(tetl+dtet)
c     *         -vi(i-1,j,k)*sin(tetl-dtet))
c . . .  Дивергенция V с котангенсом:
            dd=dd+vi(i,j,k)/(r(k)+re)*(cos(tetl)/sin(tetl))+
     *         (vi(i+1,j,k)-vi(i-1,j,k))/(r(k)+re)/dtet*0.5
            alf=dd
            fp(k)=qq
            if(alf.gt.0.)ql(k)=ql(k)+alf
            if(alf.le.0.)fp(k)=fp(k)-alf*an3(i,j,k)
 3        continue
c*******
          fp(1)=0.
          fp(lv)=0.
          do 2 k=1,np
            dx=r(k+1)-r(k)
            bbc2=(bbc(k)+bbc(k+1))/2.
            aal2=(aal(k)+aal(k+1))/2.
            aaa=(bbc2/aal2+bbc(k+1)/aal(k+1))/4.*dx
            bbb=(bbc2/aal2+bbc(k)/aal(k))/4.*dx
            if (aaa.gt.40.) then
             print 777, i,j,k,aaa
             aaa=40.
            endif
            ed(k)=exp(aaa)*aal2/dx
            if (bbb.gt.40.) then
             print 778, i,j,k,bbb
             bbb=40.
            endif
            ep(k)=exp(-bbb)*aal2/dx
            ppm=(vr(i,j,k+1)-abs(vr(i,j,k+1)))/2.
            ppp=(vr(i,j,k)+abs(vr(i,j,k)))/2.
            if((vr(i,j,k+1).lt.0.).and.(vr(i,j,k).gt.0.))ppm=0.
            ppma=(b2(k+1)-abs(b2(k+1)))/2.
            pppa=(b2(k)+abs(b2(k)))/2.
            ppmb=(b3(k+1)-abs(b3(k+1)))/2.
            pppb=(b3(k)+abs(b3(k)))/2.
            ed(k)=ed(k)-ppm-ppma-ppmb
            ep(k)=ep(k)+ppp+pppa+pppb
 2        continue
c*******
          ap(1)=0.
          bp(1)=an3(i,j,1)
          zz=ed(1)
c*******
          do 6 k=2,np
            ddx=(r(k+1)-r(k-1))/2.
            dtt=dt
            cp(k)=ddx*(1./dtt+ql(k))
            fp(k)=(an3(i,j,k)/dtt+fp(k))*ddx
            aaa=ep(k)+cp(k)+zz
            zzz=ed(k)*(cp(k)+zz)/aaa
            ap(k)=ed(k)/aaa
            bp(k)=(fp(k)+ep(k-1)*bp(k-1))/aaa
            zz=zzz
 6        continue
c*******
          hgr=r(lv)-r(np)
          ho=bk*an6(i,j,lv)/g(lv)/am3
          ppo(np)=bp(np)/(1.-ep(np)/ed(np)*ap(np))
          hod=1./ho
          ppo(lv)=ppo(np)*exp(-hgr*hod)
          an31(i,j,lv)=ppo(lv)
c*******
          do 7 jj=2,lv
            k=lv-jj+1
            ppo(k)=ap(k)*ppo(k+1)+bp(k)
            an31(i,j,k)=(ppo(k))
 7        continue
401   continue
 40   continue
  777 format(' GSO. MISTAKE  aaa ',3i4,1pe10.2)
  778 format(' GSO. MISTAKE  bbb ',3i4,1pe10.2)

      deallocate (aal,bbc,dd31,dd32,ql,
     *          fp,ap,bp,ep,ed,b2,
     *          b3,cp,dpp,pfu,ppo)


      return
      end