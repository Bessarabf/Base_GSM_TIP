      subroutine gsmo(pgl,kpars,nh,its,ids,r,i,j,dt,
     *                vr1,aal,bbc,fp,p)
      dimension pgl(kpars,nh,its,ids),r(nh),vr1(nh,its,ids),
     *          aal(nh),bbc(nh),p(nh),pg(nh),fp(nh)
      allocatable  ap(:),bp(:),ep(:),ed(:),
     *             cp(:),dpp(:),pfu(:),ppo(:)
	allocate (ap(nh),bp(nh),ep(nh),ed(nh),	
     *             cp(nh),dpp(nh),pfu(nh),ppo(nh))
	data alfs/0.10/
c     data alfs/0.25/
      np=nh-1
      do 2 k=2,np
        dk=-p(k)
        ddx=(r(k+1)-r(k-1))/2.
        ap(k)=-(ddx*dk-ddx/dt)
        bp(k)=ddx*(fp(k)+(pgl(6,k,i,j)/dt))
 2    continue
      do 6 k=2,nh
          dx=r(k)-r(k-1)
          bbc2=(bbc(k)+bbc(k-1))/2.
          aal2=(aal(k-1)+aal(k))/2.
          ep(k)=(aal2/dx-bbc2/2.)
      if(abs(ep(k)).gt.0.2)goto 10
        ed(k)=aal2/dx+bbc2/2.
          cp(k)=0.
          dpp(k)=0.
          goto 6
 10       ed(k)=0.
          cp(k)=1./ep(k)
          dpp(k)=(aal2/dx+bbc2/2.)/ep(k)
 6        continue
c         **********
          ap(1)=(r(2)-r(1))/dt/2.
          bp(1)=ap(1)*pgl(6,1,i,j)
          b00=0.
          g00=pgl(6,1,i,j)
          bp(nh)=(r(nh)-r(nh-1))/dt/2.*(pgl(6,nh,i,j))
          ap(nh)=(r(nh)-r(nh-1))/dt/2.
          fnp1=0.
      call progp(ppo,pfu,ap,bp,ep,cp,dpp,nh,b00,g00,
     *ed,fnp1)
          do 100 k=1,nh
            pgl(6,k,i,j)=(pfu(k))
            if (pgl(6,k,i,j).le.1.e-6) pgl(6,k,i,j)=1.e-3
c                    Labtam
c           vr1(k,i,j)=ppo(k)/pgl(6,k,i,j)
c                    Labtam
c                    PC
            vr1(k,i,j)=ppo(k)/pgl(6,k,i,j)
c                    PC
 100      continue
c     сглаживание
          do 1 k=2,np
            pg(k)=(1.-2.*alfs)*pgl(6,k,i,j)+
     +      alfs*(pgl(6,k-1,i,j)+pgl(6,k+1,i,j))
   1      continue
          do 3 k=2,np
            pgl(6,k,i,j)=pg(k)
   3      continue
c
c         print *,'gsmo END'
      deallocate (ap,bp,ep,ed,	
     *             cp,dpp,pfu,ppo)
          return
          end

