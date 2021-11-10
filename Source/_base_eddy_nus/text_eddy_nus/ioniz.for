      subroutine ioniz(sole,solen,rads,par,parj,gkoor,uts,dtets,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars,mass,ps)
      dimension sole(nse),solen(nse),par(kpars,nh,its),parj(nh,its),
     *       mass(30),rads(nh),gkoor(2,its,ids),ps(10)
      call ionizu(sole,solen,rads,par,gkoor,uts,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars)
      if(mass(12).eq.0)go to 1
        call ionize(par,parj,uts,dolm,dtets,delta,nh,its,kpars,ps)
    1 continue
      return
      end
