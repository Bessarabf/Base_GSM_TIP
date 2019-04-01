      subroutine qdismod(qdis,an1,an6,gkoor,g,rads,solu,nsu,delta,
     *            nh,its,ids,uts)
      dimension an1(its,ids,nh),an6(its,ids,nh),qdis(2,its,ids,nh),
     *          solu(nsu),rads(nh),gkoor(2,its,ids),g(nh)

      data pi,om/3.14159265359d0,7.272205e-5/
      cr=180./pi
      del=delta
	do i=1,its
	  do j=1,ids
          gshir=gkoor(1,i,j)/cr
          gdol=gkoor(2,i,j)/cr
          gshir=pi/2.-gshir
          f=sin(gshir)*sin(del)+cos(gshir)*cos(del)
     *      *cos(om*(uts-43200.)+gdol)
          hi=acos(f)
      
          do 1 k = 1 , nh
            temp=an6(i,j,k)
            anq=an1(i,j,k)
        
            qdis(1,i,j,k)=1.e9*dismod(anq,temp,g(k),rads(k),
     *                    solu,nsu,hi,2)
            qdis(2,i,j,k)=dismod(anq,temp,g(k),rads(k),solu,nsu,hi,1)
        
  1       continue
        end do
	end do
      return
      end
