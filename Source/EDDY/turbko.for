      subroutine turbko(eddyco,pgl,rads,kpars,n,its,ids)
      dimension eddyco(n,its),rads(n),pgl(kpars,n,its,ids)
! constants c1=cmax eddy diffusion k-t, c0=c(80km)       
      data c1,c0 /1.0e7,5.0e5/

      data s1,s2,s3/ 0.01,0.01,0.01 /
      hm=120. ! height of Kt maximum (really < 120.)
      do i=2,its-1
      do k=1,n
! zonal mean density and temperatures
         anO=0.
         anO2=0.
         anN2=0.
         tn=0.
         do j=1,ids
           anO2=anO2+pgl(1,k,i,j)
           anN2 =anN2+pgl(2,k,i,j)
           anO=anO+pgl(3,k,i,j)
           ano2=ano2+pgl(1,k,i,j)
           tn=tn+pgl(7,k,i,j)
         end do 
         sum=(anO2+anN2+anO)/ids
         Temp=tn/ids
         cMol=3.e17/sum*sqrt(Temp)
!!!!!!!!!!!
         ah=(rads(k))*1.e-5
         y=hm-ah
         y2=y*y
         if (y2.gt.5000.) then
           eddyco(k,i)=0.
         else
           if(ah.lt.hm) then
             ctd1=(c1-c0)*exp(-s2*y2)
             eddyco(k,i)=ctd1+c0*exp(-s3*y)
           else
             eddyco(k,i)=c1*exp(-s1*y2)
           end if
         end if
!!!!!!!!!!! Kt=0 if Kt<Kmol
         eddyco(k,i)=eddyco(k,i)*(1.+sign(1.,eddyco(k,i)-cMol))*0.5 
       end do
      end do
      return
      end