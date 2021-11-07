      subroutine gsoom  (an1,an11,an2,an3,vr,an6,q,
     *                  rp,r,g,lv,dt,ctd,ro,n1,n2,n,vi,vj)
c     last redaction  8.02.96
      dimension an1(n1,n2,n),an2(n1,n2,n),
     *          an3(n1,n2,n),an6(n1,n2,n),an11(n1,n2,n),
     *          vr(n1,n2,n),q(n1,n2,n),rp(n),r(n),g(n),ctd(n),
     *          ro(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n),
     *          dp(3,3)
      allocatable aal(:),bbc(:),dd12(:),dd13(:),
     *          fp(:),ap(:),bp(:),ep(:),ed(:),
     *          b2(:),b3(:),ql(:),
     *          cp(:),dpp(:),pfu(:),ppo(:),cs(:)
      allocate (aal(n),bbc(n),dd12(n),dd13(n),
     *          fp(n),ap(n),bp(n),ep(n),ed(n),
     *          b2(n),b3(n),ql(n),
     *          cp(n),dpp(n),pfu(n),ppo(n),cs(n))
      data am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
      data bk/1.38e-16/,re/6.381e08/
      data do2,dn2,do/3.6e-8,3.7e-8,2.e-8/,pi/3.141592/
      np=lv-1
      n11=n1-1
      n21=n2
      dtet=pi/(n11)
      dfi=2.*pi/n2
      gam=1.e-20
	
c*******
      do 40 i=2,n11
         tetl=dtet*(i-1)
        do 401 j=1,n21
         jp=j+1
         jl=j-1
         if(j.eq.1)jl=n21
         if(j.eq.n21)jp=1
c  *************************************************************
         cs(1)=an1(i,j,1)+an2(i,j,1)+an3(i,j,1)
	
        do 100 k=2,n
        amsv=(32.*an1(i,j,k)+an2(i,j,k)*28.+16.*an3(i,j,k))
     *       /(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))*1.67e-24
        amsn=(32.*an1(i,j,k-1)+an2(i,j,k-1)*28.+16.*an3(i,j,k-1))
     *       /(an1(i,j,k-1)+an2(i,j,k-1)+an3(i,j,k-1))*1.67e-24
            hnn=bk*an6(i,j,k-1)/g(k-1)/amsn
            hvv=bk*an6(i,j,k)/g(k)/amsv
            expp=exp(-(hnn+hvv)/(2.*hnn*hvv)*(r(k)-r(k-1)))
            cs(k)=cs(k-1)*an6(i,j,k-1)/an6(i,j,k)*expp
 100    continue
c  *************************************************************
          do 1 k=1,lv
            ams=(32.*an1(i,j,k)+an2(i,j,k)*28.+16.*an3(i,j,k))
     *         /(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))*1.67e-24
            cono2=an1(i,j,k)
            conn2=an2(i,j,k)
            cono=an3(i,j,k)
            tempp=an6(i,j,k)
            plotg=ro(i,j,k)
c  *************************************************************
          sig=0.5*(do2+dn2)
          sig1=0.5*(do2+do)
c . . . к-т диффузии O2-N2
          dp12=3./8./sig/(an1(i,j,k)+an2(i,j,k))/sig*
     *    sqrt(bk*an6(i,j,k)*(am1+am2)/am1/am2/2./pi)
c . . . к-т диффузии O2-O
          dp13=3./8./sig1/(an1(i,j,k)+an3(i,j,k))/sig1*
     *    sqrt(bk*an6(i,j,k)*(am1+am3 )/am1/am3/2./pi)
          dp(1,1)=(dp12*an2(i,j,k)+dp13*an3(i,j,k))
     *            /(an2(i,j,k)+an3(i,j,k))
c  *************************************************************
            hh=bk*tempp/g(k)
            hcp=hh/ams
            hoo=hh/am1
            hnn=hh/am2
            ho=hh/am3
            kv=k+1
            kn=k-1
            if(k.eq.lv)kv=lv
            if(k.eq.1)kn=1
            gtem=(an6(i,j,kv)-an6(i,j,kn))/(r(kv)-r(kn))
            aal(k)=ctd(k)+dp(1,1)
            bbc(k)=ctd(k)*(1./hcp+1./tempp*gtem)+
     *             dp(1,1)*(1./hoo+1./tempp*gtem)
 1        continue
c
          do 3 k=2,np
            alf=9.9e-34*exp(470./an6(i,j,k))
            bet=1.1e-34*exp(510./an6(i,j,k))
            if(q(i,j,k).ne.0.)bet=0.
            pp=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
            qq=(alf*pp+gam)*an3(i,j,k)**2
            ql(k)=q(i,j,k)+bet*an3(i,j,k)*pp
c*******
            fp(k)=qq
c . . .  Дивергенция V:
            aa=0.5/dfi/(r(k)+re)/sin(tetl)
            dd=aa*(vj(i,jp,k)-vj(i,jl,k))
            aa=0.5/dtet/(r(k)+re)/sin(tetl)
c            dd=dd+aa*(vi(i+1,j,k)*sin(tetl+dtet)
c     *         -vi(i-1,j,k)*sin(tetl-dtet))
c . . .  Дивергенция V с котангенсом:
            dd=dd+vi(i,j,k)/(r(k)+re)*(cos(tetl)/sin(tetl))+
     *         (vi(i+1,j,k)-vi(i-1,j,k))/(r(k)+re)/dtet*0.5
            alf=dd
            if(alf.gt.0.)ql(k)=ql(k)+alf
            if(alf.le.0.)fp(k)=fp(k)-alf*an1(i,j,k)
 3        continue
          fp(1)=0.
          fp(lv)=0.
          do 2 k=1,np
            bbc2=(bbc(k)+bbc(k+1))*.5
            aal2=(aal(k)+aal(k+1))*.5
            dx=r(k+1)-r(k)
            aaa=((bbc2/aal2+bbc(k+1)/aal(k+1))/4.*dx)
            bbb=((bbc2/aal2+bbc(k)/aal(k))/4.*dx)
            if(aaa.gt.40.) then
             print 777,aaa
             aaa=40.
            endif
            ed(k)=exp(aaa)*aal2/dx
            if(bbb.gt.40.) then
             print 778,bbb
             bbb=40.
            endif
            ep(k)=exp(-bbb)*aal2/dx
            ppm=(vr(i,j,k+1)-abs(vr(i,j,k+1)))/2.
            ppp=(vr(i,j,k)+abs(vr(i,j,k)))/2.
      if((vr(i,j,k+1).lt.0.).and.(vr(i,j,k).gt.0.))ppm=0.
            ed(k)=ed(k)-ppm
            ep(k)=ep(k)+ppp
 2        continue
c*******
          ap(1)=0.
          bp(1)=an1(i,j,1)
          zz=ed(1)
c*******
          do 6 k=2,np
            ddx=(r(k+1)-r(k-1))/2.
            cp(k)=ddx*(1./dt+ql(k))
            fp(k)=(fp(k)+an1(i,j,k)/dt)*ddx
            aaa=ep(k)+cp(k)+zz
            zzz=ed(k)*(cp(k)+zz)/aaa
            ap(k)=ed(k)/aaa
            bp(k)=(fp(k)+ep(k-1)*bp(k-1))/aaa
            zz=zzz
 6        continue
c*******
          hgr=r(lv)-r(np)
          hoo=bk*an6(i,j,lv)/g(lv)/am1
          ppo(np)=bp(np)/(1.-ap(np)*(ep(np)/ed(np)))
          hood=1./hoo
          ppo(lv)=ppo(np)*exp(-hgr*hood)
          an11(i,j,lv)=ppo(lv)
c*******
          do 7 jj=2,lv
            k=lv-jj+1
            ppo(k)=ap(k)*ppo(k+1)+bp(k)
            an11(i,j,k)=(ppo(k))
 7        continue
401   continue
 40   continue
      deallocate (aal ,bbc ,dd12 ,dd13 ,
     *          fp ,ap ,bp ,ep ,ed ,
     *          b2 ,b3 ,ql,
     *          cp ,dpp ,pfu ,ppo ,cs )
  777 format(' GSOO MISTAKE aaa=',1pe10.2)
  778 format(' GSOO MISTAKE bbb=',1pe10.2)
      return
      end