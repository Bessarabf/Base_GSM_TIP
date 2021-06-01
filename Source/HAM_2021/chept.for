      function chept(x,hi)
c     . . . CHEPMEN's function by Titheridge
      data pi/3.141592/
	
c     definition f-formula
      ch(a,b,hi)=a*(1./cos(b*hi)-0.834)
	if(hi.gt.pi) hi=pi
      hig=hi/pi*180.
      a=3.88/x**1.143
      sq=sqrt(x)
      b=1.0123-1.454/sq
      if(hig.lt.89.) then
       d=ch(a,b,hi)
       chept=1./cos(hi-d)
       return
      else
       z=sin(hi)*x
	
       sqz=sqrt(z)
       pok=x-z
       if(pok.gt.30.)pok=30.
       ex=exp(pok)
       prom=2.507*sqz+.903/sqz
       hip=pi-hi
       d=ch(a,b,hip)
       chept=ex*prom-1./cos(hip-d)
      endif
      return
      end

