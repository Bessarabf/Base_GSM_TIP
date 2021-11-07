! version with massive ionization from HAMMONIA 22/03/2019
      subroutine ioniz_AEham(sole,solen,rads,par,parj,gkoor,uts,dtets,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars,mass,ps,E0,FAE,j)
!!!!!
      USE mo_ham_gsm
!!!!!     
      dimension sole(nse),solen(nse),par(kpars,nh,its),parj(nh,its),
     *       mass(30),rads(nh),gkoor(2,its,ids),ps(10)
      dimension E0(its,*),FAE(its,*)
!!!! new mass for compatibility with HAMMONIA particles ionization
     *         ,peion(4,nh0,its0)

!!!   coefficient for precipitation ionization (like in ionizv.for)
      dimension cc(4)
!!      data cc /1.,1.,0.6,1./
! NO partical ionization = 0
      data cc/1.09,1.,0.61,0./

      call ionizu(sole,solen,rads,par,gkoor,uts,
     *       dolm,ddolgs,delta,nh,its,ids,nse,kpars)
      if(mass(12).eq.1) then
!  our model precipitatios
         call ionize(par,parj,uts,dolm,dtets,delta,nh,its,kpars,ps)
       else if(mass(12).ge.2.and.mass(12).le.5) then
!  our model+Zhang&Paxton or Vorobjov&Yagodkina auroral precipitation
         call ionize_AE_ham(par,parj,peion,uts,dolm,ddolgs,dtets,
     *                  delta,nh,its,kpars,ps,E0,FAE)
         do kp =1,4 
           do i=1,its
             do k=1,nh
               par(kp+12,k,i)=par(kp+12,k,i)+peion(kp,k,i)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               parj(k,i)=peion(1,k,i)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
             end do
           end do
        end do
         
       else if(mass(12).eq.6) then 
!  HAMMONIA ionization before mlev (120 km) + Z&P above
         call ionizE_AE_ham(par,parj,peion,uts,dolm,ddolgs,dtets,
     *                  delta,nh,its,kpars,ps,E0,FAE)	
!         do i= 1,its  
!          write(*,*)'ionization',i,j,zpiongsm(:,i,j)*1.e-6
!          write(*,*)'ionization orig',i,j,peion(1,:,i)
!         enddo
         do kp=1,4 
           do i=1,its
             do k=1,mlev+5
               par(kp+12,k,i)=par(kp+12,k,i)+
     *                        cc(kp)*zpiongsm(k,i,j)*1.e-6
             end do
! smoothiing profile ionization rate
             par(kp+12,mlev+6,i)=par(kp+12,mlev+6,i)+
     *                          (peion(kp,mlev+6,i)+cc(kp)*
     *                           zpiongsm(mlev+6,i,j)*1.e-6)*0.5
             do k=mlev+7,nh
               par(kp+12,k,i)=par(kp+12,k,i)+peion(kp,k,i)
             end do
           end do
        end do
      else
	   print*,' without electron precipitation!'

      end if
      return
      end
	  