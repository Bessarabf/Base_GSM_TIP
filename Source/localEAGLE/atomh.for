      subroutine atomh(pgl,kpars,nh,its,ids,gkoor,rads,
     *                 fs,fa,ap,uts,day)
c     . . . расчет атомарного водорода по MSIS до 520 км
      integer day
      dimension pgl(kpars,nh,its,ids),tm(2),gkoor(2,its,ids)
     *         ,apm(7),rads(nh)
!      dimension dm(8) ! MSIS 90
      dimension dm(9)  ! MSIS 2000
      data om/7.27e-5/,pi/3.14159/
      iyd=80*1000+day
c     vys=rads(nh)/1.e5
      do l=1,7
        apm(l)=ap
      end do
c     write(*,*) fs,fa,ap
      do 88 i = 1 , its
        do 88 j = 1 , ids
           fig=gkoor(1,i,j)
           fig=90.-fig
           dgeo=gkoor(2,i,j)/180.*pi
           dol=gkoor(2,i,j)
           tault=uts+dgeo/om
           tault=tault/3600.
  111      if(tault.gt.24.) then
            tault=tault-24.
            go to 111
           end if
           tau2=uts
           if(tau2.gt.86400.) tau2=tau2-86400.
           do k=15,nh
             vys=rads(k)*1.e-5
             ! MSIS2000
             call gtd7(iyd,tau2,vys,fig,dol,tault,fs,fa,apm,
     *                  1,dm,tm)
             pgl(5,k,i,j)=dm(7)
           end do
88       continue
      close(25)
      return
      end
