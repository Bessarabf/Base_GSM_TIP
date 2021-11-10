      subroutine timol(pgl1,kpars,nh,its,ids,dts,vim,vir,vid)
      real vim(nh,its,ids),vid(nh,its,ids),vir(nh,its,ids),nu0,
     *     pgl1(kpars,nh,its,ids)
      data bk/1.38e-16/,nu0/0.9e-9/,ami/30./,ae/1.6e-24/
  900 format(' ',10g12.3)
      a0=2*7.69e-19/3/bk
      c0=ami*ae*nu0/6/bk
      igp=its-1
      do 3 j = 1 , ids
       do 1 ig=2,its-1
        do 2 i=1,nh
      if(pgl1(9,i,ig,j).le.0.) print22,pgl1(9,i,ig,j),i,ig,j
  22   format(' error , te=',e10.3,' i=',i4,' ig=',2i4)
          dvr=pgl1(10,i,ig,j)-vir(i,ig,j)
          dvt=pgl1(11,i,ig,j)-vim(i,ig,j)
          dvf=pgl1(12,i,ig,j)-vid(i,ig,j)
          dvr2=dvr*dvr
          dvt2=dvt*dvt
	

          dvf2=dvf*dvf
          a=a0*pgl1(6,i,ig,j)/(pgl1(9,i,ig,j)**1.5)
          sn=pgl1(1,i,ig,j)+pgl1(2,i,ig,j)+0.45*pgl1(3,i,ig,j)
          b=nu0*sn/2.
          c=c0*sn*(dvr2+dvt2+dvf2)
          r1=dts*(a*pgl1(9,i,ig,j)+b*pgl1(7,i,ig,j)+c)
          r2=1.+dts*(a+b)
          pgl1(8,i,ig,j)=(pgl1(8,i,ig,j)+r1)/r2
          if(pgl1(8,i,ig,j).le.pgl1(7,i,ig,j)) pgl1(8,i,ig,j)=
     *   pgl1(7,i,ig,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          pgl1(8,i,ig,j)=pgl1(7,i,ig,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	     
    2   continue
    1 continue
    3   continue
      return
      end
                                                                                                                                                                                                                                                                           
