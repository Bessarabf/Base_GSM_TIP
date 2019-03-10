      subroutine aprete(cqm,fa,hsl,jz,nol,ns,nu,nvn,pga,pun,qs3,
     *                  tsl,wm,qp)
      integer jz,nol,ns,nu(*),nvn
      real cqm,fa,hsl(nvn,*),pga,pun,qs3(*),tsl(nvn,*),wm(*),qp
      real aqfn(23),aqfs(23),qfn,qfs
      call cqfns(cqm,fa,hsl,jz,nol,nu,nvn,qs3,tsl,qfn,qfs)
      call allqfns(ns,qfn,qfs,wm,aqfn,aqfs)
      call ctoen(aqfn,aqfs,ns,pga,pun,wm,qp)
      return
      end