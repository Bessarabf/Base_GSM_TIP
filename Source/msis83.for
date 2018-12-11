      subroutine msis83(day,hl,alat,along,f,fbar,key,n,ng,z,zg,
     *t,doxy,do2,dn2,dh,tg,dgoxy,dgo2,dgn2)
      dimension doxy(n),dn2(n),dh(n),do2(n),t(n),z(n),
     *dgoxy(n),dgn2(n),dgo2(n),tg(n),zg(n)
      common/bldela/ap1,ap2,ap3,ap4,ap12av,ap36av,utsec
      dzeta(z1,za1)=(z1-za1)*(6.35677e3+za1)/(6.35677e3+z1)
      xzdzet(z1,dzeta0)=-dzeta(z1,1.165e2)/dzeta0+1.e0
      dif(tl1,tz1,c1,c2,z1)=(tl1/tz1)**c1*exp(c2*z1)
      dif1(c3,c4,t01,tb1,tc1,td1,xz1)=
     *c3*exp(c4*(t01*(xz1-1.e0)+3.33333e-1*tb1*(xz1**3-1.e0)+
     *2.e-1*tc1*(xz1**5-1.e0)+1.4286e-1*td1*(xz1**7-1.e0)))
c
      cosd1=cos(1.721e-2*(day+4.485e0))
      cosday=cos(1.7214e-2*(day+7.692e0))
      sinlt=sin(2.618e-1*hl)
      sinlt2=sin(5.236e-1*hl)
      sinlt3=sin(7.854e-1*hl)
      coslt=cos(2.618e-1*hl)
      coslt2=cos(5.236e-1*hl)
      coslt3=cos(7.854e-1*hl)
      delf=f-fbar
      delfav=fbar-1.50e2
c                Labtam
c         pust0=0.e0
c         pust1=1.e0
c         pust2=1.e0
c         pust3=1.e0
c         zet1=z(j1)
c         zet2=z(j1)
c           zetg1=zg(j1)
c           zetg2=zg(j1)
c                Labtam
c
      x=sin(alat)
      x2=x*x
      x3=x2*x
      x4=x2*x2
      x5=x4*x
      x6=x5*x
      x12=sqrt(1.e0-x2)
      x32=x12**3
      p10=x
      p11=x12
      p20=1.5e0*x2-5.e-1
      p21=3.e0*x*x12
      p22=-3.e0*x2+3.e0
      p30=2.5e0*x3-1.5e0*x
      p31=x12*(7.5e0*x2-1.5e0)
      p32=1.5e1*(x-x3)
      p33=1.5e1*x32
      p40=4.375e0*x4-3.75e0*x2+3.75e-1
      p41=x12*(1.75e1*x3-7.5e0*x)
      p42=-5.25e1*x4+6.0e1*x2-7.5e0
      p43=1.05e2*x*x32
      p50=7.875e0*x5-8.75e0*x3+1.875e0*x
      p51=x12*(3.9375e1*x4-2.625e1*x2+1.875e0)
      p52=-1.575e2*x5+2.1e2*x3-5.25e1*x
      p60=1.44375e1*x6-1.96875e1*x4+6.5625e0*x2-3.125e-1
      p61=x12*(8.6625e1*x5-7.875e1*x3+1.3125e1*x)
c 50
      p63=x32*(1.7325e3*x3-4.725e2*x)
c
      tlgrad=1.466e+1*(1.e0+p20*(1.377e-2+2.039e-3*delfav)+9.448e-2*p40
     *+1.872e-3*delfav+coslt2*(2.851e-2*p22+7.593e-3*p42-5.954e-3*cosday
     **p32)+sinlt2*(2.718e-3*p22+1.922e-3*p42+6.611e-3*cosday*p32))
      texos=1.035e3*(1.e0+2.394e-2*p20-1.616e-3*p40+1.631e-2*p60+delfav*
     *(3.067e-3-6.771e-6*delfav-4.669e-4*p20)+delf*(1.778e-3-5.460e-6*
     *delf)+
     * 1.098e-2*cos(1.721e-2*(day-3.762e+1))+1.424e-2*cos(3.442e-2*
     *(day-1.304e+2))-(1.561e-1*p10+3.303e-2*p30)*(1.e0+3.041e-3*delfav
     *+delf*(1.778e-3-5.460e-6*delf))*cosday-1.853e-2*p10*cos(3.443e-2*
     *(day-6.667e0))+(1.e0+5.83e-3*delfav+delf*(1.778e-3-5.46e-6*delf))*
     *((-1.12e-1*p11-6.432e-3*p31+1.146e-2*p21*cosday)*coslt+(-1.273e-1*
     *p11+1.191e-4*p31-7.107e-3*p21*cosday)*sinlt+(-2.736e-3*p22-
     *1.175e-3*p42+3.948e-3*p32*cosday)*coslt2+(9.537e-3*p22-1.087e-3*
     *p42+1.897e-3*p32*cosday)*sinlt2+(8.287e-4*p33+(-1.034e-4*p43+
     *1.232e-4*p63)*cosday)*coslt3+(1.15e-3*p33-(2.196e-4*p43+3.038e-6*
     *p63)*cosday)*sinlt2)+dela(4.031e-5,1.723e-1,2.379e-1)*(7.634e-3+
     *9.267e-3*p20+2.703e-3*p40+(1.518e-3*p21+1.728e-3*p41+6.895e-4*p61)
     **(1.e0-3.833e-1*p10*cosday)*cos(along+6.961e1)+(2.262e-4*p10+
     *2.225e-3*p30-6.283e-3*p60)*cos(7.272e-5*( utsec-1.221e+4)))+(1.e0+
     *5.521e-3*delfav)*((3.078e-3*p21+2.418e-3*p41+5.23e-3*p61)*
     *cos(along)+(-2.636e-2*p21-2.261e-3*p41+3.664e-3*p61)*sin(along))-
     *(1.277e-2*p10+4.166e-2*p30+2.886e-2*p50)*(1.e0+5.264e-3*delfav)*
     *(1.e0-3.914e-1*p10)*cos(7.272e-5*(utsec-2.940e+4))+(5.751e-4*p11+
c 75
     *3.794e-4*p31+9.714e-4*p51)*cos(7.272e-5*utsec-along))
      texos=texos+1.035e3*((-4.516e-3*p11-2.752e-3*p31+1.116e-3*p51)*
     *sin(7.272e-5*utsec-along)+(5.116e-4*p32+5.452e-4*p52)*
     *cos(7.272e-5*(utsec-9.608e+2)+2.e0*along))
      tl=3.807e+2*(1.e0+7.269e-4*delfav+1.229e-2*cos(3.443e-2*(day-
     *7.147e+1))-3.53e-2*p10*cosday+p11*(-1.106e-2*coslt+1.148e-2*sinlt)
     *+(-8.919e-3*p22-2.63e-3*p42+3.152e-4*p52*cosday)*coslt2
     *+(1.275e-2*p22-6.404e-4*p42+2.593e-3*p52*cosday)*sinlt2
     *+dela(4.031e-5,1.723e-1,2.379e-1)*(8.026e-3+6.997e-3*p20))
      t0=1.806e+2*(1.e0-7.325e-2*p20+(1.276e-1*p10+8.189e-2*p30)*cosday
     *+(9.927e-3*p22-5.105e-4*p42)*coslt2
     *+(1.578e-3*p22+2.647e-3*p42)*sinlt2)
      z0=9.455e+1*(1.e0-2.637e-2*p20+1.924e-4*delfav+4.722e-2*p10*
     *  cosday
     *+(6.48e-3*p22-2.568e-3*p42+1.314e-3*p52*cosday)*coslt2
     *+(4.491e-3*p22-4.101e-4*p42+1.543e-3*p52*cosday)*sinlt2)
      tr=1.591e-1*(1.0e0+1.006e0*p20-5.88e-1*p10*cosday
     *+(-6.160e-2*p22+3.401e-2*p42)*coslt2
     *-(1.082e-1*p22+4.640e-2*p42)*sinlt2)
c
      denol=8.603e10*exp(-3.312e-2*p10-1.486e-1*p20-9.209e-2*p40
     *+delfav*(2.694e-3-8.728e-6*delfav)+5.235e-4*delf+
     *dela(2.443e-5,1.088e-1,3.413e-1)*(-1.412e-2-6.598e-2*p20-
     *  1.338e-2*
     *p40+(-4.474e-3*p21-5.619e-3*p41-4.022e-3*p61)*cos(along+5.488e+1)
     *+(-1.104e-2*p10+1.07e-2*p30)*cos(7.272e-5*(utsec+1.067e4)))+
     *8.679e-2*cos(1.721e-2*(day+1.763e1))+
     *1.348e-1*cos(3.443e-2*(day-1.059e2))+
     *(3.722e-1*p10+1.017e-2*p30)*(1.e0+5.235e-4*delf)*cosd1+
     *(1.e0+1.381e-3*delfav+5.235e-4*delf)*
     *((-6.652e-2*p11-4.13e-3*p31+2.637e-2*p21*cosd1)*coslt+
     *(6.518e-2*p11-3.165e-2*p31+5.197e-2*p21*cosd1)*sinlt+
     *(1.111e-2*p22+3.585e-3*p42)*coslt2-
     *(2.866e-2*p22+9.286e-4*p42)*sinlt2+
     *p33*(8.635e-4*coslt3-1.186e-3*sinlt3))+
     *(1.e0+2.77e-3*delfav)*((-3.803e-3*p21-4.561e-3*p41-1.025e-2*p61)*
     *cos(along)+(7.889e-2*p21+9.617e-3*p41-1.42e-2*p61)*sin(along))+
     *(7.36e-3*p10-1.68e-1*p30-1.337e-1*p50)*(1.e0+5.056e-3*delfav)*
     *(1.e0-3.833e-1*p10)*cos(7.272e-5*(utsec+1.383e+4))-
     *(3.885e-3*p32+2.254e-3*p52)*cos(7.272e-5*(utsec+4.747e3)+2.e0*
     *along))
      deno2l=3.165e10*exp(1.999e-1*p20-1.683e-3*delfav+
     *dela(+3.766e-5,1.723e-1,2.379e-1)*(8.551e-3+2.474e-2*p20)+
     *5.307e-2*cos(1.721e-2*(day-1.516e1))+
     *2.809e-2*cos(3.443e-2*(day-8.337e1))+6.394e-2*p10*cosday+
     *p11*(-3.207e-2*coslt+4.736e-2*sinlt)-
     *p22*(1.508e-2*coslt2+4.651e-2*sinlt2)-
     *(1.828e-2*p21+8.966e-3*p41)*cos(along)-
     *(6.266e-2*p21+4.113e-2*p41)*sin(along))
      denn2l=3.296e11*exp(1.670e-3*delfav+
     *7.413e-2*cos(1.721e-2*(day+2.211e1))+
     *3.789e-2*cos(3.443e-2*(day-1.351e2))-9.98e-2*p10*cosday-
     *p11*(1.576e-2*coslt+2.674e-2*sinlt)+
     +(1.851e-2*p22+6.376e-3*p42)*coslt2+
     *(-3.798e-2*p22+3.335e-3*p42)*sinlt2+
     *p33*(1.406e-4*coslt3-1.832e-3*sinlt3))
      denhl=2.185e+5*exp(-4.768e-2*p10-2.074e-1*p20-1.614e-1*p40-
     *1.368e-2*delfav+delf*(-6.316e-3+4.094e-5*delf)+
     *dela(2.904e-5,1.723e-1,2.379e-1)*(-2.834e-2-1.29e-2*p20+
     *2.03e-2*p40+(7.843e-3*p10+6.565e-3*p30)*
     *cos(7.272e-5*(utsec-3.165e3)))+
     *5.53e-2*cos(1.721e-2*(day-1.425e2))+1.854e-2*cos(3.443e-2*(day-
     *7.791e+1))+(4.272e-1*p10+1.18e-1*p30)*(1.e0+delf*(-6.316e-3+
     *4.094e-5*delf))*cos(1.721e-2*(day+8.994e0))+
     *(1.e0+delf*(-6.316e-3+4.094e-3*delf)-6.697e-3*delfav)*
     *(p11*(2.101e-1*coslt+2.78e-1*sinlt)+
     *p31*(2.937e-2*coslt+2.191e-2*sinlt)+
     *p21*cos(1.721e-2*(day+8.994e0))*(2.573e-2*coslt+3.728e-2*sinlt)+
     *p22*(1.72e-2*coslt2-1.111e-2*sinlt2)-
     *p33*(1.366e-3*coslt3+2.501e-4*sinlt3))-
     *(5.769e-3*p21+1.823e-2*p41+1.236e-2*p61)*cos(along)+
     *(5.762e-2*p21+5.666e-3*p41-1.707e-3*p61)*sin(along)+
     *(-5.067e-2*p10+9.676e-2*p30+1.232e-1*p50)*
     *cos(7.272e-5*(utsec+1.791e+4)))
      dzeta0=dzeta(z0,1.165e+2)
      tdelt=texos-tl
      sigma=tlgrad/tdelt
      ta=texos-tdelt*exp(3.5019e0*sigma)
      tagrad=(texos-ta)*sigma*1.00108e0
      t12=t0+tr*(ta-t0)
      t0=1.e0/t0
      tb=1.e0/ta-t0
      tc=dzeta0*tagrad/(ta*ta)
      td=6.6666e-1*tc-3.11111*tb+7.11111*(1.e0/t12-t0)
      tc=5.e-1*tc-tb-2.e0*td
      tb=tb-tc-td
      c1h=1.1362e0/(texos*sigma)
      c1n2=c1h*2.8e1
      c1o2=c1h*3.2e1
      c1o=c1h*1.6e1
      c10=c1h*2.895e1
      c2h=-1.13623e0/texos
      c2n2=c2h*2.8e1
      c2o2=c2h*3.2e1
      c2o =c2h*1.6e1
      c20 =c2h*2.895e1
      c3n2=dif(tl,ta,c1n2,c2n2,-3.5019e0)
      c3o2=dif(tl,ta,c1o2,c2o2,-3.5019e0)
      c3o =dif(tl,ta,c1o ,c2o ,-3.5019e0)
      c3h =dif(tl,ta,c1h ,c2h ,-3.5019e0)
      c30 =dif(tl,ta,c10 ,c20 ,-3.5019e0)
      c4h=1.1375e0*dzeta0
      c4n2=c4h*2.8e1
      c4o2=c4h*3.2e1
      c4o=c4h*1.6e1
      c40=c4h*2.895e1
      xz95=xzdzet(9.5e1,dzeta0)
      xz100=xzdzet(1.00e2,dzeta0)
      xz105=xzdzet(1.05e2,dzeta0)
      eddyav=dif1(c30,c40,t0,tb,tc,td,xz105)
      eddyn2=dif1(c3n2,c4n2,t0,tb,tc,td,xz105)/eddyav*denn2l
      eddyo2=dif1(c3o2,c4o2,t0,tb,tc,td,xz105)/eddyav*deno2l
      eddyo=dif1(c3o,c4o,t0,tb,tc,td,xz100)/eddyav*denol
      eddyh=dif1(c3h,c4h,t0,tb,tc,td,xz95)/eddyav*denhl
      rato2=eddyn2/eddyo2
      rato =eddyn2/eddyo
      rath =eddyn2/eddyh
      assign 1 to lable1
      assign 6 to lable2
      assign 9 to lable3
      if (key.gt.1) assign 4 to lable3
      if (key.eq.3) assign 5 to lable2
c
        do 9 j=1,n
          j1=n-j+1
          go to lable1,(1,2)
   1      if (z(j1).ge.1.165e2) assign 2 to lable1
          if (z(j1).ge.1.165e2) go to 2
            xz=xzdzet(z(j1),dzeta0)
            tz=1.e0/(t0+tb*(xz**2)+tc*(xz**4)+td*(xz**6))
            t(j)=tz
            difn2=dif1(c3n2,c4n2,t0,tb,tc,td,xz)
            difo2=dif1(c3o2,c4o2,t0,tb,tc,td,xz)
            difo=dif1(c3o,c4o,t0,tb,tc,td,xz)
            difh=dif1(c3h,c4h,t0,tb,tc,td,xz)
            eddy0=dif1(c30,c40,t0,tb,tc,td,xz)
            go to 3
  2       dz=dzeta(z(j1),1.20e2)
            tz=texos-tdelt*exp(-sigma*dz)
            t(j)=tz
            difn2=dif(tl,tz,c1n2,c2n2,dz)
            difo2=dif(tl,tz,c1o2,c2o2,dz)
            difo =dif(tl,tz,c1o,c2o,dz)
            difh =dif(tl,tz,c1h,c2h,dz)
            eddy0=dif(tl,tz,c10,c20,dz)
   3      dn2(j)=den(denn2l,difn2,tl,tz,0.e0,2.8e1,eddyn2,eddy0,pust1,
     *    pust2,zet1,zet2,pust3,pust0,z(j1),1.e0)
          do2(j)=den(deno2l,difo2,tl,tz,pust0,3.2e1,eddyo2,eddy0,
     *    2.683e-1,rato2,zet1,1.152e2,6.376e0,0.e0,z(j1),1.e0)
          doxy(j)=den(denol,difo,tl,tz,0.e0,1.6e1,eddyo,eddy0,8.218e-2,
     *    rato,z(j1),1.326e2,4.026e1,-3.346e1,7.166e1,8.161e0)
          dh(j)=den(denhl,difh,tl,tz,-4.e-1,1.e0,eddyh,eddy0,2.329e-6,
     *    rath,z(j1),1.223e2,2.386e1,-1.700e1,7.593e1,2.020e0)
          go to lable3,(4,9)
  4       go to lable2,(5,6)
c
  5       if (j1.gt.ng) go to 9
   6      if (zg(j1).ge.1.165e2) assign 7 to  lable2
          if (zg(j1).ge.1.165e2) go to 7
            xz=xzdzet(zg(j1),dzeta0)
            tz=1.e0/(t0+tb*(xz**2)+tc*(xz**4)+td*(xz**6))
            tg(j)=tz
            difn2=dif1(c3n2,c4n2,t0,tb,tc,td,xz)
            difo2=dif1(c3o2,c4o2,t0,tb,tc,td,xz)
            difo=dif1(c3o,c4o,t0,tb,tc,td,xz)
            difh=dif1(c3h,c4h,t0,tb,tc,td,xz)
            eddy0=dif1(c30,c40,t0,tb,tc,td,xz)
            go to 8
  7       dz=dzeta(zg(j1),1.20e2)
            tz=texos-tdelt*exp(-sigma*dz)
            tg(j)=tz
            difn2=dif(tl,tz,c1n2,c2n2,dz)
            difo2=dif(tl,tz,c1o2,c2o2,dz)
            difo =dif(tl,tz,c1o,c2o,dz)
            difh =dif(tl,tz,c1h,c2h,dz)
            eddy0=dif(tl,tz,c10,c20,dz)
   8      dgn2(j)=den(denn2l,difn2,tl,tz,pust0,2.8e1,eddyn2,eddy0,pust1,
     *    pust2,zetg1,zetg2,pust3,0.e0,zg(j1),1.e0)
          dgo2(j)=den(deno2l,difo2,tl,tz,pust0,3.2e1,eddyo2,eddy0,
     *    2.68e-1,rato2,zetg1,1.152e2,6.376e0,0.e0,zg(j1),1.e0)
          dgoxy(j)=den(denol,difo,tl,tz,0.e0,1.6e1,eddyo,eddy0,8.218e-2
     *    , rato,zg(j1),1.326e2,4.026e1,-3.346e1,7.166e1,8.161e0)
  9   continue
      return
      end

