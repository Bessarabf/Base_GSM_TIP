      subroutine sfera(alt,tzbaz,cn1,cn,pltn,tn)
      dimension cn(6),cn1(6)
      dimension cnzr(6)
      dimension alf(6),amk(6),am(6),alfa(6)
      data alf/0.,0.,0.,0.,0.17,-0.38/
      data amk/28.,32.,0.,16.,4.,40./
      data zbaz,zr,zm,zs/125.,100.,85.,45./
      data tm,ts/185.,270./
      data cnzr/1.e19,2.e18,0.,1.e18,1.e17,7.e13/
      data zp,hp/95.,2.5/
      data pi/3.1415927/,sq/1.7320508/
      data r,gzer,re/8.31434e3,9.80665,6371.02/
c
      recomb(x,y)=(1.+exp(-(y-zp)/hp))/(1.+exp(-(x-zp)/hp))
      cuba(x)=x**3+1
      g(x)=gzer*(re/(re+x))**2.
      cubic(x)=atan((2.*x-1.)/sq)/sq-alog(1.-3.*x/(1.+x)**2.)/6.
c
      deltaz=zbaz-zm
      a=(tzbaz/tm-1.)**(1./3.)
      do 16 i=1,6
      cnzr(i)=cn1(i)*1.e6
   16 continue
      st4=deltaz*1.e3/(a*r*tm)
      st3=(zm-zs)*2.*1.e3/(r*sqrt(tm*ts)*pi)
      st2=sqrt(ts/tm)
      tnn=tzbaz
      if(alt.ge.85.) go to 2
            do 12 i=1,6
               am(i)=28.86
               alfa(i)=0.
   12       continue
            x=0.
            tn=((tm+ts)+(tm-ts)*cos(pi*(alt-zm)/(zs-zm)))/2.
            do  3 i=1,6
              cub=cubic(a) -cubic(x)
              cub=cub*g((zm+zbaz)/2.)
              cnpr=cnzr(i)*(tnn/tm)**(1.+alfa(i))*exp
     *        (am(i)*st4*cub)
              if(i.eq.4) cnpr=cnpr*recomb(zm,zr)
             cub=st3*am(i)*g((alt+zm)/2.)
              cn(i)=cnpr*tm/tn
              if(abs(alt-zs).lt.1.e-4) go to 14
                  atg=atan(st2*tan(pi*(zm-alt)/2./(zm-zs)))
                  if(alt.lt.zs)atg=atg+pi
                  go to 15
   14             atg=pi/2.
   15         continue
              cn(i)=cn(i)*exp(cub*atg)
                  if(i.eq.4) cn(i)=cn(i)*recomb(alt,zm)
    3       continue
            go to 13
    2       if(alt.ge.125.) go to 1
              do 8 i=1,6
                 am(i)=28.86
                  alfa(i)=0.
    8         continue
                  x=(alt-zm)/deltaz*a
                  tn=tm*cuba(x)
                  do 10 i=1,6
                 cub=cubic(a) -cubic(x)
              cn(i)=cnzr(i)*(tnn/tn)**(1.+alfa(i))
              cn(i)=cn(i)*exp(am(i)*st4*g((alt+zbaz)/2.)*cub)
                    if(i.eq.4) cn(i)=cn(i)*recomb(alt,zr)
   10             continue
    1         continue
   13 continue
      do 17 i=1,6
        cn(i)=cn(i)/1.e6
   17 continue
      pltn=0.
      do 5 i=1,6
           pltn=pltn+cn(i)*amk(i)*1.66e-24
    5 continue
      return
      end

