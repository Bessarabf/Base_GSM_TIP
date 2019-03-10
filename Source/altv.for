cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     N2 VIBRATIONAL QUANTA APPROXIMATION
c     Pavlov, A.V., The role of vibrationally excited nitrogen in the
c     formation of the mid-latitude negative ionospheric storms,
c     Annales Geophysicae, 1994, Vol.12, P. 554-564.
c
c     EQUATION (B.3)
c
c     INPUT:
c     h(nn)  - ALTITUDE(KM)
c
c     h(1)<140 km !!!  h(nn)>400 km !!!!  h(n)-h(n-1)<20 km !!!
c
c     yo(nn) - O NUMBER DENSITY(CM-3)
c     tn(nn) - TEMPERATURE
c     te(nn) - ELECTRON TEMPERATURE
c     cne(nn)- ELECTRON NUMBER DENSITY (CM-3)
c     nn - the number of altitude points
c
c    !!!!!  nn<200    !!!!!
c
c     OUTPUT:
c     alfa(nn) - VIBRATIONAL QUANTA
c     tv(nn)   - VIBRATIONAL TEMPERATURE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine altv(h,yo,tn,te,cne,alfa,tv,nn)
      dimension h(nn),yo(nn),tn(nn),te(nn),cne(nn),alfa(nn),
     *tv(nn),atv1(5),btv1(5),ctv1(5),atv2(5),btv2(5),ctv2(5)
      dimension clam(23),dl(23),w1(23),w2(23)
      double precision w(200),x1,w0,dn2,tdif,tvt,hn2,cl,cl1,w00
        data atv1/2.8,-2.745,-3.073,-4.,-4.469/
        data btv1/1.221e-3,4.454e-3,4.596e-3,4.994e-3,5.166e-3/
        data ctv1/-8.613e-8,-5.592e-7,-5.771e-7,-6.288e-7,-6.497e-7/
        data atv2/4.159,3.817,3.696,3.379,3.155/
        data btv2/7.042e-4,8.051e-4,8.303e-4,8.864e-4,9.232e-4/
        data ctv2/-4.166e-8,-4.814e-8,-4.959e-4,-5.301e-8,-5.505e-8/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     atv1,btv1,ctv1 for Te<4000 K,  atv2,btv2,ctv2 for Te>4000 K
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc  THE CALCULATION OF THE FREQUENCY AND LAMDA  cccc
        do 1 n=1,nn
        te0=te(n)
        w0=0.
        if(te0.lt.1200.) go to 4
        if(te0.gt.4000.) go to 3
        do 2 m=1,5
        x1=atv1(m)+btv1(m)*te0+ctv1(m)*te0*te0-16.d0
        x1=10.d0**x1
2       w0=w0+x1
        go to 4
3     continue
        do 5 m=1,5
        x1=atv2(m)+btv2(m)*te0+ctv2(m)*te0*te0-16.d0
        x1=10.d0**x1
5       w0=w0+x1
4     continue
c
       w(n)=w0*cne(n)
c
        g=981./(1.+h(n)/6378.)**2
        t=tn(n)
        hn2=826.e5*t/g/28.
        x1=yo(n)/1.d9
        dn2=9.69d7*t**0.724/x1
        tdif=hn2*hn2/dn2
        tvt=0.107*exp(-69.9/t**0.33)*x1
        tvt=1.d0/tvt
        clam(n)=2.d0*dsqrt(tdif/tvt)
        dl(n)=dsqrt(tdif*tvt)
1     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccc
         w0=0.
         w00=0.
         do 6 n=2,nn
         j=nn+1-n
         j1=j+1
        cl=clam(j)
        cl1=clam(j1)
        x1=(cl-cl1)/2.d0
        cl=dexp(-cl)
        cl1=dexp(-cl1)
        w0=w0+(w(j)*cl+w(j1)*cl1)*x1
        w00=w00+(w(j)/cl+w(j1)/cl1)*x1
         w1(j)=w0
         w2(j)=w00
6     continue
         w1(nn)=0.
         w2(nn)=0.
ccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 7 n=1,nn
        x1=-3353./tn(n)
        teta=dexp(x1)
        cl=clam(n)
        cl1=dexp(cl)
        x1=1.d0/cl1
        alfa(n)=teta+dl(n)/cl*(w0*(cl1-x1)-cl1*w1(n)+x1*w2(n))
        cl=alfa(n)
	if(cl.le.0) then
	   print*,'altv n=',n,'cl=', cl
	   cl=10.d0
	end if 
        tv(n)=-3353./dlog(cl/(1.d0+cl))
7     continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
