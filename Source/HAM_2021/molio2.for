      subroutine molio2 (pgl,pgi,rads,vir,kpars,nh,its,dts,
     *     ids,ddolgs,ins,dtet,ntr)
c
c     расчет ион. m+ c переносом по вертикали
c

      dimension pgl (kpars,nh,its,ids),rads(nh),vir(nh,its,ids),
     *     pgi(ins,nh,its,ids)
     
      allocatable dk(:),qs(:),r1(:),r2(:),
     *     p(:),ce(:),co(:)
	

      real nu0,l1,l2,lamo
      data nu0/1.e-9/,ami/30./,ae/1.6e-24/,el/1.6e-20/,b0/0.3/,
     *     pi/3.141592/,re/6371.02e5/,bk/1.38e-16/,alfa/4.2e-7/,
     *     l1/1.6e-11/,l2/6.e-13/,ge/980.665/
  771 format(' molion2: kpars,nh,its,ids,ddolgs,dtet,ntr :'/
     *       ' ',4i6,2g12.4,i4)

      allocate (dk(nh),qs(nh),r1(nh),r2(nh),
     *     p(nh),ce(nh),co(nh))

	cr=180./pi
      ame=ae*ami/el
      gkm=ae*ami*ge/bk
      nmd=360./ddolgs+0.1
c     nd=dolm/ddolgs+1
      itd=its-1
      do 1 j=1,ids
         do 2 i=2,itd
c
c        расчет коэф. зависящих от широты
c
         tet=dtet*(i-1)/cr
         dip=1.
         if(tet.eq.0.)dip=pi/2.
         if(tet.eq.pi)dip=-pi/2.
         if(tet.eq.pi/2.)dip=0.
         if(dip.ne.1.)go to 3
           t=tan(tet)
           dip=atan(2./t)
    3 continue
         si=sin(dip)
         si2=si*si
         ct=cos(tet)
         bb=b0*sqrt(1.+3*ct*ct)
         be=ame/bb
         bnu=be*nu0
c
c        расчет коэф. диффузии и фотохимии
c
         do 4 k=1,nh
c       **************
           pr=(pgl(1,k,i,j)+pgl(2,k,i,j))/2.+pgl(3,k,i,j)/3.
           g=bnu*pr
           g2=g*g
           rj=ae*ami*nu0*pr
           dk(k)=bk*(si2+g2)/(1.+g2)/rj
           ce(k)=pgl(6,k,i,j)
           qs(k)=pgl(13,k,i,j)+pgl(14,k,i,j)+pgl(15,k,i,j)
c          if(k.ge.ntr+2) go to 5
           if(k.gt.ntr+2) go to 5
             alt=alfa*(300./pgl(9,k,i,j))
             co(k)=pgl(16,k,i,j)/(l1*pgl(1,k,i,j)+l2*pgl(2,k,i,j))
             ce(k)=pgl(6,k,i,j)+co(k)
             qs(k)=qs(k)+pgl(16,k,i,j)
             p(k)=alt*ce(k)
             go to 6
    5      continue
           alt=alfa*(300./pgi(5,k,i,j))
           ce(k)=ce(k)+pgi(1,k,i,j)
           vrn=pgl(10,k,i,j)
           vtn=pgl(11,k,i,j)
           vln=pgl(12,k,i,j)
           vior=pgi(2,k,i,j)
           viot=pgi(3,k,i,j)
           viol=pgi(4,k,i,j)
           tn=pgl(7,k,i,j)
c
           tim=tn
           if(k.ge.16)tim=pgi(6,k,i,j)
           l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
           l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
c          te=pgl(9,k,i,j)
c          l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c          l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c          l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol)
c          l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol)
c
           qs(k)=qs(k)+(l1*pgl(1,k,i,j)+l2*pgl(2,k,i,j))*
     *     pgi(1,k,i,j)
           p(k)=alt*ce(k)
c       if(p(k).lt.ry)print776,k,alt,ce(k),pgi(5,k,i,j)
c 776 format(' molion2: k,alt,ce(k),pgi(5,k,i,j)'/' ',4g13.4)
    6      continue
    4   continue
c
c       расчет коэфф. при производных
c
        n1=nh-1
        do 7 k=2,n1
        h2=rads(k+1)-rads(k-1)
        cf1=(pgl(9,k+1,i,j)+pgl(8,k+1,i,j)-pgl(9,k-1,i,j)-
     *  pgl(8,k-1,i,j))/h2
        cf2=(pgi(1,k+1,i,j)-pgi(1,k-1,i,j))/h2
        cf2=cf2*pgi(5,k,i,j)/ce(k)
        r1(k)=dk(k)*(pgl(8,k,i,j)+pgl(6,k,i,j)*pgl(9,k,i,j)/
     *  ce(k))
        r2(k)=dk(k)*(cf1+cf2+gkm)-vir(k,i,j)
    7  continue
c
c      расчет коэф. на границах
c
        h2=rads(3)-rads(1)
        r1(1)=dk(1)*(pgl(8,1,i,j)+pgl(6,1,i,j)*pgl(9,1,i,j)/ce(1))
c       cf1=pgl(9,3,i,j)+pgl(8,3,i,j)-4*(pgl(9,2,i,j)+pgl(8,2,i,j))
c       cf1=cf1+3*(pgl(9,1,i,j)+pgl(8,1,i,j))/h2
c       cf2=(pgi(1,3,i,j)-4*pgi(1,2,i,j)+3*pgi(1,1,i,j))/h2
c       cf2=cf2*pgi(5,1,i,j)/ce(1)
c       r2(1)=dk(1)*(cf1+cf2+gkm)-vir(1,i,j)
        r2(1)=0.
c       h2=rads(nh)-rads(nh-2)
        h1=rads(nh)-rads(nh-1)
        r1(nh)=dk(nh)*(pgl(8,nh,i,j)+pgl(6,nh,i,j)*pgl(9,nh,i,j)/
     *  ce(nh))
c       cf1=(pgl(9,nh,i,j)+pgl(8,nh,i,j))*3+pgl(9,nh-2,i,j)+
c    *  pgl(8,nh-2,i,j)-4*(pgl(9,nh-1,i,j)+pgl(8,nh-1,i,j))
c       cf1=cf1/h2
c       cf2=(pgi(1,nh-2,i,j)-4*pgi(1,nh-1,i,j)+3*pgi(1,nh,i,j))/h2
c       cf2=cf2*pgi(5,nh,i,j)/ce(nh)
        cf1=(pgl(9,nh,i,j)+pgl(8,nh,i,j))
     *  -(pgl(9,nh-1,i,j)+pgl(8,nh-1,i,j))
        cf1=cf1/h1
        cf2=(pgi(1,nh,i,j)-pgi(1,nh-1,i,j))/h1
        cf2=cf2*pgi(5,nh,i,j)/ce(nh)
        r2(nh)=dk(nh)*(cf1+cf2+gkm)-vir(nh,i,j)

c
c       граничное условие для концентрации
c
        pgl(6,1,i,j)=(pgl(6,1,i,j)+qs(1)*dts)/(1.+p(1)*dts)
c
c       прогонка
c
        call gsmo (pgl,kpars,nh,its,ids,rads,i,j,dts,
     *  vir,r1,r2,qs,p)
    2   continue
    1   continue
        print *,'molion2 END'
	deallocate (dk,qs,r1,r2,
     *     p,ce,co)
        return
        end

