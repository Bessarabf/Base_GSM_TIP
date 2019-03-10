      subroutine ioniz_AE(sole,solen,rads,par,parj,gkoor,uts,dtets,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars,mass,ps,E0,FAE)
      dimension sole(nse),solen(nse),par(kpars,nh,its),parj(nh,its),
     *       mass(30),rads(nh),gkoor(2,its,ids),ps(10)
      dimension E0(its,*),FAE(its,*)
      call ionizu(sole,solen,rads,par,gkoor,uts,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars)
      if(mass(12).eq.1) then
	!  our model
         call ionize(par,parj,uts,dolm,dtets,delta,nh,its,kpars,ps)
       else if(mass(12).ge.2.and.mass(12).le.5) then
      !  our model+Zhang&Paxton or Vorobjov&Yagodkina auroral precipitation
         call ionize_AE(par,parj,uts,dolm,ddolgs,dtets,
     *                   delta,nh,its,kpars,ps,E0,FAE)
	else
	   print*,' without electron precipitation!'

      end if
      return
      end
	  