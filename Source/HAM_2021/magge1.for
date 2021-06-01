      subroutine g11t21(lpaz,nj,y,jww,ih,ind,ux)
      dimension y(nj)
      jww=11021
      ind=1
      ux=0.
      lp=iabs(lpaz)
      r=float(lp)
      ih=j52t62(nj,nj,r,y)
      s=y(ih)
      lt=ifix(s)
      u=y(ih+1)
      lu=ifix(u)
      if(lpaz.lt.0) then
        la=lp-lt
        lb=lu-lp
        if(la.gt.lb) then
          ih=ih+1
          lpaz=lu
          goto 3
        end if
        lpaz=lt
        goto 3
      end if
        ind=2
        ux=(r-s)/(u-s)
    3 continue
      return
      end
      function g11t22(gir,ifs,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
      dimension gir(lr),pars(kp,nhs,nvs,lr),tes(nvs)
      p=pars(ifs,i,j,l)
      if(kp.eq.8.and.ifs.le.3) goto 1
        g=alog10(p)
        goto 2
    1   g=0.43429*p
    2 continue
      g11t22=g
      return
      end
      function g11t23(gir,ifs,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
      dimension gir(lr),pars(kp,nhs,nvs,lr),tes(nvs)
      g11t23=pars(ifs,i,j,l)
      return
      end
      function g11t24(gir,ifs,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
      dimension gir(lr),pars(kp,nhs,nvs,lr),tes(nvs)
      p=pars(ifs,i,j,l)*.01
      g11t24=p
      return
      end
      function g11t25(gir,ifs,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
      dimension gir(lr),pars(kp,nhs,nvs,lr),tes(nvs)
      if(lgmgr.eq.0) then
        p=g11t24(gir,ifs,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
      else
        t=tes(j)
        f=gir(l)
        r=g11t31(f,t)
        s=sin(r)
        c=cos(r)
        ik=11
        st=g11t24(gir,ik,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
        ik=12
        sf=g11t24(gir,ik,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
        if(ifs.eq.11) then
          p=st*c-sf*s
        else
          p=sf*c+st*s
        end if
      end if
      g11t25=p
      return
      end
      subroutine g11t26(g5s,gir,ifs,kpars,lchose,lgmgr,lpaz,lr,ng,
     *                  nh,nhh,nim,nv,pars,tes,x,jww,psdv)
      dimension gir(lr),pars(kpars,nh,nv,lr),tes(nv),x(ng),
     *          psdv(nim,nim)
      jww=11026
      g=float(iabs(lpaz))
      goto (1,2),lchose
    1 continue
        tg=g
        goto 3
    2 continue
        dg=g
    3 continue
      ks=0
      do 8 ig=1,ng
        goto (4,5),lchose
    4   continue
          dg=x(ig)
          goto 6
    5   continue
          tg=x(ig)
          if(ig.eq.1) then
            dg=30.
          else if(ig.eq.ng) then
            dg=30.
          else
            dg=g
          end if
    6   continue
        call g52t61(ks,dg,tg,dm,tm)
        i=j52t62(lr,lr,dm,gir)
        j=j52t62(nv,nv,tm,tes)
        k=i+1
        l=j+1
        r=gir(i)
        cx=(dm-r)/(gir(k)-r)
        r=tes(j)
        cy=(tm-r)/(tes(l)-r)
        r=1.-cy
        s=1.-cx
        a=r*s
        b=r*cx
        c=cy*s
        d=cy*cx
        do 7 ih=1,nhh
          pa=g5s(gir,ifs,ih,j,kpars,i,lgmgr,lr,nh,nv,pars,tes)
          pb=g5s(gir,ifs,ih,j,kpars,k,lgmgr,lr,nh,nv,pars,tes)
          pc=g5s(gir,ifs,ih,l,kpars,i,lgmgr,lr,nh,nv,pars,tes)
          pd=g5s(gir,ifs,ih,l,kpars,k,lgmgr,lr,nh,nv,pars,tes)
          psdv(ig,ih)=a*pa+b*pb+c*pc+d*pd
    7   continue
    8 continue
      return
      end
      function g11t27(gir,ifs,i,j,kp,l,lgmgr,lr,nhs,nvs,pars,tes)
      dimension gir(lr),pars(kp,nhs,nvs,lr),tes(nvs),am(4)
      data am/5.32e-23,4.67e-23,2.66e-23,1.67e-24/
      s=0.
      do 1 k=1,4
        s=s+pars(k,i,j,l)*am(k)
    1 continue
      g11t27=alog10(s)
      return
      end
      function g11t31(f,t)
      double precision a,b,c,e,ps,r,s,u,v
      data s/1.976573d-1/,c/9.802712d-1/
      data ps/1.74532925199432d-2/
      e=dble(t)*ps
      u=dble(f)*ps
      a=dcos(e)
      b=dcos(u)
      r=a*b*s
      a=dsin(e)
      v=c*a
      e=v+r
      a=dsin(u)
      r=s*a
      b=datan2(r,e)
      g11t31=sngl(b)
      return
      end

