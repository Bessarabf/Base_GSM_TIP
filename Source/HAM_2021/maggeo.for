      subroutine g523(gir,hs,ifs,kpars,lchose,lgmgr,lpaz,lr,
     *                nhds,nh,nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
      external g11t22,g11t23,g11t24,g11t25,g11t27
      dimension gir(lr),hs(nh),pars(kpars,nh,nv,lr),tes(nv),
     *          psdv(nim,nim),x(nim),y(nim)
      jww=523
      goto (1,1,1,1,1,1,2,2,2,3,4,4,1,1,1,1,5),ifs
    1 continue
        call g52t11(g11t22,gir,hs,ifs,kpars,lchose,lgmgr,lpaz,lr,
     *              nhds,nh,nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
        goto 7
    2 continue
        call g52t11(g11t23,gir,hs,ifs,kpars,lchose,lgmgr,lpaz,lr,
     *              nhds,nh,nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
        goto 7
    3 continue
        call g52t11(g11t24,gir,hs,ifs,kpars,lchose,lgmgr,lpaz,lr,
     *              nhds,nh,nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
        goto 7
    4 continue
        call g52t11(g11t25,gir,hs,ifs,kpars,lchose,lgmgr,lpaz,lr,
     *              nhds,nh,nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
        goto 7
    5 continue
        call g52t11(g11t27,gir,hs,ifs,kpars,lchose,lgmgr,lpaz,lr,
     *              nhds,nh,nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
    7 continue
      return
      end
c
      subroutine g52t11(g5s,gir,hs,ifs,kpars,lchose,lgmgr,lpaz,lr,
     *                  nhh,nh,nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
       dimension gir(lr),hs(nh),pars(kpars,nh,nv,lr),tes(nv),
     *           psdv(nim,nim),x(nim),y(nim)
      external g5s
       jww=52011
       call g52311(g5s,gir,hs,ifs,kpars,lgmgr,lpaz,lr,nh,nhh,nim,nv,
     *               pars,tes,jww,ni,nj,psdv,x,y)
       return
       end

      subroutine g52311(g5s,gir,hs,ifs,kpars,lgmgr,lpaz,lr,nh,nhh,
     *                  nim,nv,pars,tes,jww,ni,nj,psdv,x,y)
      dimension gir(lr),hs(nh),pars(kpars,nh,nv,lr),tes(nv),
     *          psdv(nim,nim),x(nim),y(nim)
      jww=52311
      ni=lr+1
      nj=nv
      do 1 j=1,nj
        y(j)=tes(j)
    1 continue
      call g11t21(lpaz,nh,hs,jww,ih,ind,ux)
      do 5 l=1,lr
        x(l)=gir(l)
        do 4 j=1,nj
          fa=g5s(gir,ifs,ih,j,kpars,l,lgmgr,lr,nh,nv,pars,tes)
          if(ind.eq.2) then
            ig=ih+1
            fb=g5s(gir,ifs,ig,j,kpars,l,lgmgr,lr,nh,nv,pars,tes)
            fa=fa+(fb-fa)*ux
          end if
          psdv(l,j)=fa
    4   continue
    5 continue
      r=x(2)-x(1)
      x(ni)=x(lr)+r
      do 6 j=1,nj
        psdv(ni,j)=psdv(1,j)
    6 continue
      if(lgmgr.eq.1) then
         call g52t41(ni,nim,nj,x,y,jww,psdv)
      else
         call g52041(ni,nim,nj,x,y,jww,psdv)
      end if
      return
      end

