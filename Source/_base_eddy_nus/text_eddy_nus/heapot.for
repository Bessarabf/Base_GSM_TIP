      subroutine heapot(pgl,pgi,an1,an2,an3,an6,anco2,vr,vi,vj,
     *                  vim,vid,vir,ctd,solu,nsu,
     *                  rads,g,gkoor,kpars,ins,nh,its,ids,
     *                  delta,uts,mass)
c
c     . . . calculate and write
c     sourse heating in labtam.file
c     1- uv
c     2- euv & precipitation
c     3- IR-cooling
c     4- Joul
c     5- convection transport
c     6- compress - decompress
c     7- condactivity cooling
c     8- viscosity heating
c     9-
c     ne & nd - number 'edenic i desyatkow w mass
c     nd may be .le. ne, ne .le. 8
c     nd is first, ne is last writing parameters.
c     then you can use sub quad
c
c
     
      character*4 lass(1:9),lamp(1:4)
      dimension an1(its,ids,nh),an2(its,ids,nh),an3(its,ids,nh),
     *          an6(its,ids,nh),vr(its,ids,nh),vi(its,ids,nh),
     *          vj(its,ids,nh),vim(nh,its,ids),vid(nh,its,ids),
     *          vir(nh,its,ids),pgl(kpars,nh,its,ids),
     *          pgi(ins,nh,its,ids),solu(nsu),
     *          gkoor(2,its,ids),rads(nh),g(nh),
     *          ctd(nh),anco2(its,ids,nh)
      allocatable h(:,:,:),an(:,:,:),alt(:)
     *         
      data lass/' uv ','euv ',' ik ',' qj ','dkon','pdv ','dtep',
     *          ' dnu',' rcv'/,re/6.371e 8/,bk/1.38e-16/
     *     lamp/'conh','pdvv','pdvh','conv'/
c
      idsp=ids+1
      nht=26
      allocate(h(its,ids,nh),an(its,ids+1,nht))
      
c
      sq=-1.
      sq1=-1.
      mab=iabs(mass)
      nd=mab/10
      if(nd.eq.0) nd=mab*11
      ne=mab-nd*10
      if(ne.lt.0.or.ne.gt.8) then
       print 100,mass
       return
      end if
c     . . . write dencity
      do 71 i=1,its
       do 72 j=1,ids
        do 73 k=1,nh
         rc=2.5*(an1(i,j,k)+an2(i,j,k))+
     *      1.5*an3(i,j,k)
         h(i,j,k)=rc*bk
   73   continue
   72  continue
   71 continue
      call rewrit(an,h,nh,its,ids,nht)
      call write8(an,lass(9),uts,nht,its,idsp,sq,sq1)
c     . . . non adiabatic source calculating
      if(nd.le.4) then
       nstop=4
       if(ne.lt.4) nstop=ne
        do 3 i=nd,nstop
         mis=isign(i,mass)
         call joupot(h,pgl,pgi,ctd,rads,g,solu,gkoor,
     *               kpars,ins,nh,its,ids,nsu,uts,delta,
     *               an1,an2,an3,an6,anco2,vr,vi,vj,vim,vid,vir,mis)
c        call zinpac(h,rads,g,gkoor,nh,its,ids,
c    *       ctd,vim,vid,vir,an1,an2,an3,an6,vi,vj,vr,mis)
c        . . .  line only with zinpac!
c        lass(i)=lamp(i)
         call bongl(h,nh,its,ids)
         call rewrit(an,h,nh,its,ids,nht)
c     . . .
      sq=-1.
      sq1=-1.
         if(.not.( mass.gt.0.)) go to 23
          do 24 k=1,nht
           alt(k)=re+rads(k)
   24     continue
c         call quadro(sq,an,alt,26,its,ids)
          call quad  (sq,sq1,an,alt,nht,its,ids,1,nht)
   23    continue
         call write8(an,lass(i),uts,nht,its,idsp,sq,sq1)
    3   continue
       end if
c     . . . dynamic sourse and condactiviti heating
      sq=-1.
      if(ne.lt.5) go to 90
       nbeg=5
       if(nd.gt.5) nbeg=5
       do 4 i=nbeg,ne
        mis=isign(i,mass)
        call dinpac(h,rads,g,gkoor,nh,its,ids,
     *              ctd,vim,vid,vir,an1,an2,
     *              an3,an6,vi,vj,vr,mis)
        call bongl(h,nh,its,ids)
        call rewrit(an,h,nh,its,ids,nht)
         call write8(an,lass(i),uts,nht,its,idsp,sq,sq1)
    4  continue
   90 continue
      nraz=ne-nd+1
      print 101,nraz
  100 format(' sub heapac'/
     *' parameter mass(18)=',i4,' uncorrect!',
     *' pacage turn off')
  101 format(' sub heapac  writing',i3,' parameters')
      deallocate(h,an)
      return
      end

c
c     sub. calculate heat sourse for Gsm tip
c     It's one of fore sourse
c     1 - convection transport
c     2 - compress
c     3 - condactivity heating
c     4 - viscosity heating
      subroutine dinpac(h,rads,g,gkoor,nh,its,ids,
     *                  ctd,vim,vid,vir,
     *                  an1,an2,an3,an6,vi,vj,vr,m)
      dimension an1(its,ids,nh),
     *          an2(its,ids,nh),an3(its,ids,nh),
     *          an6(its,ids,nh),h(its,ids,nh),
     *          vj(its,ids,nh),vi(its,ids,nh),
     *          vr(its,ids,nh),rads(nh)
     *,         ctd(nh),g(nh),gkoor(2,its,ids)
     *,         vim(nh,its,ids),vid(nh,its,ids)
     *,         vir(nh,its,ids),alyam(3),ro(3)
      data re/6.371e8/,pi/3.1415926/,bk/1.38e-16/,
     *    am1,am2,am3/53.12e-24,46.51e-24,26.56e-24/
       mab=iabs(m)-4
      n0=nh-1
      itsm1=its-1
      dteta=pi/itsm1
      ddol=2.*pi/ids
c     . . .   latitude and longitude in geom. coor. frame
      do 1 i=2,itsm1
       teta=dteta*(i-1)
       sin t=sin(teta)
       sin p=sin(teta+dteta)
       sin m=sin(teta-dteta)
      do 2 j=1,ids
      jm=j-1
      jp=j+1
      if(j.eq.1) jm=ids
      if(j.eq.ids) jp=1
      do 3 k=2,n0
c     .  .  . condactivity and dencity
        do 20 k1=1,3
         kv=k-2+k1
          al=18.6*an6(i,j,kv)**0.84*an1(i,j,kv)
          al=al+27.2*an6(i,j,kv)**0.8*an2(i,j,kv)
          al=al+67.1*an6(i,j,kv)**0.71*an3(i,j,kv)
          alyam(k1)=al/(an1(i,j,kv)+an2(i,j,kv)+
     *              an3(i,j,kv))
          rocp=(3.5*(an1(i,j,kv)+an2(i,j,kv))+2.5*
     *               an3(i,j,kv))*bk
          alyam(k1)=alyam(k1)+rocp*ctd(kv)
c
          ro(k1)=an1(i,j,kv)*am1+
     *           an2(i,j,kv)*am2+
     *           an3(i,j,kv)*am3
   20   continue
        sum=(an1(i,j,k)+an2(i,j,k)+an3(i,j,k))*bk
        rc=(2.5*(an1(i,j,k)+an2(i,j,k))+1.5*an3(i,j,k))*bk
c       . . . step on altitude
        drad=rads(k+1)-rads(k)
        drad1=rads(k)-rads(k-1)
c     . . . altitude  rk
        rk=rads(k)+re
        go to(10,11,12,13),mab
   14   return
c         heat  convection transport
   10   con1=vr(i,j,k)-abs(vr(i,j,k))*.5
        con2=vr(i,j,k)+abs(vr(i,j,k))*.5
        tvk=con1*(an6(i,j,k+1)-an6(i,j,k))/drad
        tvk=tvk+con2*(an6(i,j,k)-an6(i,j,k-1))/drad1
        con1=vi(i,j,k)-abs(vi(i,j,k))*.5
        con2=vi(i,j,k)+abs(vi(i,j,k))*.5
        tvi=con1*(an6(i+1,j,k)-an6(i,j,k))/(rk*dteta)
        tvi=tvi+con2*(an6(i,j,k)-an6(i-1,j,k))/(rk*dteta)
        con1=vj(i,j,k)-abs(vj(i,j,k))*.5
        con2=vj(i,j,k)+abs(vj(i,j,k))*.5
        tvj=con1*(an6(i,jp,k)-an6(i,j,k))/(rk*ddol*sin t)
        tvj=tvj+con2*(an6(i,j,k)-an6(i,jm,k))/(rk*ddol*sin t)
        h(i,j,k)=-(tvk+tvi+tvj)*86400.
        go to 3
c       . . . compress
   11   davl0=(vr(i,j,k+1)-vr(i,j,k-1))/(drad1+drad)
        davl1=(vi(i+1,j,k)-
     *         vi(i-1,j,k))/dteta*.5+vi(i,j,k)*cos(teta)/sin t
        davl2=(vj(i,jp,k)-vj(i,jm,k))/ddol*.5
        davl=davl0+(davl1+davl2)/(rk*sin t)
        h(i,j,k)=-davl*sum*an6(i,j,k)/rc*86400.
        go to 3
c        condactivity cooling
   12   tep=(alyam(3)+alyam(2))*(an6(i,j,k+1)-an6(i,j,k))/drad
        tep=(tep-(alyam(2)+alyam(1))*(an6(i,j,k)-an6(i,j,k-1))/drad1)/
     *             (drad+drad1)
        h(i,j,k)=tep/rc*86400.
        go to 3
c        viscosity heating
   13   anu=3.34e-6*an6(i,j,k)**0.71
        vqi=(vi(i,j,k+1)-vi(i,j,k))/drad
        vqj=(vj(i,j,k+1)-vj(i,j,k))/drad
        vq=(vqi*vqi+vqj*vqj)*anu
        h(i,j,k)=vq/rc*86400.
    3  continue
      h(i,j,nh)=h(i,j,n0)
      h(i,j,1)=h(i,j,2)
    2 continue
    1 continue
  223 format(2x,1p8e9.2)
      return
      end
c
      subroutine joupot(q,pgl,pgi,ctd,rads,g,solu,
     *                  gkoor,kpars,ins,nh,its,ids,nsu,uts,
     *                  del,an1,an2,an3,an6,anco2,vr,vi,vj,
     *                  vim,vid,vir,m)
c    . . . sub calculate heart source that as:
c    solar uv and euv heating, IR- cooling, joule and etc.
c
c  cold  flux, djoulp. koh. co2                         12.10.92 BFS
c
      dimension q(its,ids,nh),g(nh),pgl(kpars,nh,its,ids),
     *          rads(nh),solu(nsu),ctd(nh),sni(6),gkoor(2,its,ids)
     *         ,an1(its,ids,nh),an2(its,ids,nh),an3(its,ids,nh)
     *         ,an6(its,ids,nh),anco2(its,ids,nh),vr(its,ids,nh)
     *         ,vi(its,ids,nh),vj(its,ids,nh)
     *         ,pgi(ins,nh,its,ids)
     *         ,vim(nh,its,ids),vid(nh,its,ids)
     *         ,vir(nh,its,ids)
c
      data pi/3.1415926/,om/7.272205e-5/,bk/1.38e-16/
c
c   . . . эффективность для зимы
     *     r0,r00/3.48,2.94/
c   . . . эффективность=0.6
c    *     r1,r2/3.48,2.94/,
c   . . . эффективность=0.4
c    *     r1,r2/2.32,1.96/
c   ..... efficency = 0.45
     *     r1,r2/2.61,2.205/
      mab=iabs(m)
      n0=nh-1
      ddol=2.*pi/ids
      dteta=pi/(its-1)
      do 21 i=1,its
       do 22 j=1,ids
c             geograf longitude and latitude
        cr=pi/180.
        g shir=gkoor(1,i,j)*cr
        g dol=gkoor(2,i,j)*cr
        g shir=pi/2.-g shir
c             zenit angle
        f=sin(g shir)*sin(del)+cos(gshir)*cos(del)*
     *    cos(om*(uts-43200.)+g dol)
        hi=acos(f)
        key=1
        q(i,j,1)=0.
        do 1 k=2,n0
         ano=an1(i,j,k)
         and=an2(i,j,k)
         antr=an3(i,j,k)
         con no=pgl(4,k,i,j)
         tem=pgl(7,k,i,j)
         viam=pgi(3,k,i,j)
         viad=pgi(4,k,i,j)
         viar=pgi(2,k,i,j)
         vimm=vim(k,i,j)
         vimd=vid(k,i,j)
         vimr=vir(k,i,j)
         vin=vi(i,j,k)
         vjn=vj(i,j,k)
         vrn=vr(i,j,k)
         cona=pgi(1,k,i,j)
         conm=pgl(6,k,i,j)
         ti=pgl(8,k,i,j)
         te=pgl(9,k,i,j)
         conco2=anco2(i,j,k)
         go to(10,11,12,13), mab
   14    return
   10  continue

         qi=dis mod(ano,tem,g(k),rads(k),solu,nsu,hi,key)*ano*0.3
         go to 2
   11    r=r1
         if(k.gt.24) r=r2
         qi=r*(pgl(13,k,i,j)+pgl(14,k,i,j)+
     *        pgl(15,k,i,j)+pgl(16,k,i,j))*1.e-11
         go to 2
   12    qi=cold(ano,and,antr,tem,conco2,rads(k),g(k),ctd(k))
         cik53=anoik(con no,antr,tem)
         qi=qi+cik53
         go to 2
c*
   13    pol=(tem+pgi(6,k,i,j))*.5
         call freac(pol,cona,conm,sni)
c . . . Дополнительный джоуль от флуктуаций Эл. поля kiss=1
         kiss=0
         call djoulp(qi,vimm,vimd,vimr,viam,viad,viar,
     *               vin,vjn,vrn,ano,and,antr,sni,kiss)
         qi=qi*2.
    2    continue
c        edinicy izmerenia?
         if(.not.(m.lt.0)) go to 3
           rc=(2.5*(ano+and)+1.5*antr)*bk
           q(i,j,k)=qi/rc*86400.
           go to 1
    3    continue
c        q(i,j,k)=alog10(qi)
         q(i,j,k)=qi
    1   continue
        q(i,j,nh)=q(i,j,n0)
        q(i,j,1)=q(i,j,2)
   22  continue
   21 continue
      return
      end
c
      subroutine rewrit(anew,an,nh,its,ids,nht)
      dimension an(its,ids,nh),anew(its,ids+1,nht)
      do 1 i=1,its
       do 2 k=1,nht
         do 3 j=1,ids
          anew(i,j,k)=an(i,j,k)
    3    continue
         anew(i,ids+1,k)=an(i,1,k)
    2  continue
    1 continue
      return
      end
c
      subroutine write8(an,ias,uts,n,n1,n2,sq,sq1)
      character*4 ias
      dimension an(n1,n2,n)
c     . . . writing three-dimensional array GSM
c                                      29.01.91
      time=uts/3600.
      write(8    ) ias,time,sq
c 100 format(a4,f8.3,1pe10.3)
      write(8    )an
c 101 format(1p100(e8.1))
      print 102,ias,sq,sq1
  102 format(' write8 ',a4,' sq=',1pe9.2,'sq1=',1pe9.2)
      return
      end
c
      subroutine quad(sq,sq1,q,r,n,n1,n2,ki,ke)
c     thri-dimensional quadrature method
c     ki - initial point
c     ke - finish point
c
      dimension q(n1,n2,n),r(n),dim(50,100),sr(30),
     *          si(50),sj(100)
     *         ,col(50),fi(100)
      data pi/3.14159/
      dtet=pi/(n1-1)
      dfi=pi/(n2-1)*2.
      ne=(n1+1)/2
      n21=n2-1
      do 11 i=1,n1
       col(i)=dtet*(i-1)
   11 continue
      do 12 j=1,n2
       fi(j)=dfi*(j-1)
   12 continue
c    . . .   norse semisphere
      do 1 i=1,ne
       do 2 j=1,n2
        do 3 k=ki,ke
         sr(k)=r(k)**2*q(i,j,k)
    3   continue
        call trap1(sr,r,n,sq,ki,ke)
        dim(i,j)=sq
    2  continue
    1 continue
      do 4 j=1,n2
       do 5 i=1,ne
        tet=dtet*(i-1)
        si(i)=dim(i,j)*sin(tet)
    5  continue
       call trap(si,col,ne,sq)
       sj(j)=sq
    4  continue
       call trap(sj,fi,n2,sq)
c    . . .   souse semisphere
      do 21 i=ne,n1
       do 22 j=1,n2
        do 23 k=ki,ke
         sr(k)=r(k)**2*q(i,j,k)
   23   continue
        call trap1(sr,r,n,sq1,ki,ke)
        dim(i,j)=sq1
   22  continue
   21 continue
      do 24 j=1,n2
       do 25 i=ne,n1
        tet=dtet*(i-1)
        si(i)=dim(i,j)*sin(tet)
   25  continue
       call trap(si,col,ne,sq1)
       sj(j)=sq1
   24  continue
       call trap(sj,fi,n2,sq1)
       return
       end