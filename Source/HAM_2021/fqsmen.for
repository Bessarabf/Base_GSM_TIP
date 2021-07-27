      subroutine fqsmen(an1,an2,an3,an6,nh,its,i,j,solu,q,
     *                  gkoor,rads,delta,ids,nsu,uts)
      dimension an1(its,ids,nh),q(its,ids,nh),
     *          solu(nsu),an2(its,ids,nh),
     *          rads(nh),gkoor(2,its,ids),
     *          an3(its,ids,nh),an6(its,ids,nh)
      data pi,om/3.14159265359d0,7.272205e-5/,g/981./
      cr=180./pi
      del=delta
      gshir=gkoor(1,i,j)/cr
      gdol=gkoor(2,i,j)/cr
      gshir=pi/2.-gshir
      f=sin(gshir)*sin(del)+cos(gshir)*cos(del)
     * * cos(om*(uts-43200.)+gdol)
      hi=acos(f)
      key=2
      do 1 k = 1 , nh
        pk=an1(i,j,k)+an2(i,j,k)+an3(i,j,k)
        temp=an6(i,j,k)
        anq=an1(i,j,k)
        
         q(i,j,k)=1.e9*dismod(anq,temp,g,rads(k),solu,nsu,hi,key)
  1    continue
      return
      end

