      subroutine degap(nvar,nc,nl,dtet,dfi,sfti,psf,pst,gamma,
     *           delta,psi,u,nui,rads,nh,na)
      dimension rads(nh),gamma(nl,nc),delta(nl,nc),psi(nl,nc),
     *          pst(nl,nc),psf(nl,nc),sfti(nl,nc),u(nui)
      data pi/3.14159265359/,re/6371.02e5/
      df=(dfi+dfi)/180.*pi
      dt=(dtet+dtet)/180.*pi
      df=1./df
      dt=1./dt
      r=re+rads(na)
      s1=(re+re)/r
c      s1=re/r
      k=nc-1
      k1=(nc+1)/2
      k2=k1-3
      k3=k1+3
      do5i=2,k
        if(i.eq.k1)goto5
          ip=i+1
          im=i-1
          if(i.le.k2)tet=im*dtet/180.*pi
          if(i.ge.k3)tet=(i-5)*dtet/180.*pi
          if(i.gt.k2.and.i.lt.k3)goto1
            s=s1*sin(tet)*cos(tet)
c            ct=cos(tet)
c            s=s1*(1.+3.*ct*ct)*sin(tet)*.5/ct
            goto2
    1     continue
            if(i.lt.k1)l=i-1
            if(i.lt.k1)goto2
              l=nc-i
    2       continue
            do4j=1,nl
              jp=j+1
              jm=j-1
              if(j.eq.1)jm=nl
              if(j.eq.nl)jp=1
              gamma(j,i)=(sfti(jp,i)-sfti(jm,i))*df
              psi(j,i)=0.
              if(nvar.ne.2)psi(j,i)=(psf(jp,i)-psf(jm,i))*df
              if(i.gt.k2.and.i.lt.k3)goto3
                delta(j,i)=(sfti(j,ip)-sfti(j,im))*dt/s
                if(nvar.ne.2)psi(j,i)=psi(j,i)-(pst(j,ip)-pst(j,im))*
     *          dt/s
                goto4
    3         continue
              du=u(l+1)-u(l-1)
              if(i.gt.k1)du=-du
              delta(j,i)=(sfti(j,ip)-sfti(j,im))/du
              if(nvar.ne.2)psi(j,i)=psi(j,i)-(pst(j,ip)-pst(j,im))/du
    4       continue
    5 continue
c     do6i=1,nl
c       do6j=5,k2
c         l=nc-j+1
c         gamma(i,j)=gamma(i,j)+gamma(i,l)
c         gamma(i,l)=gamma(i,j)
c         delta(i,j)=delta(i,j)+delta(i,l)
c         delta(i,l)=delta(i,j)
c         psi(i,j)=psi(i,j)+psi(i,l)
c         psi(i,l)=psi(i,j)
c   6 continue
      return
      end
