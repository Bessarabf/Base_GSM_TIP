      subroutine hplm(plm,hms,h,int)
      dimension plm(int),am(7)
      data am   /32.,28.,16.,30.,16.,1.,4./,re/6371.02e5/
      data alf/.8/
c     data alf/.7/
c     data alf/.5/
cc    data alf/.1/
c     data alf/.2/
c     data alf/0./
c     data alf/.95/
      bolc=1.38041e-16
      at1=1.66e-24
      g0=980.665
      ge=g0/(1+h/re)**2
      tkg=bolc*plm(8)/ge/at1
      rh=h-hms
      hr=(re+h)/(re+hms)
      i=1
    1 if(i.gt.7) go to 2
        hi=tkg/am(i)
        argum=rh/hi*hr
        if(argum.ge.100.)argum=100.
        if(i.eq.6)then
          ho=tkg/am(3)
          ho=rh/ho*hr
          plm(i)=plm(i)*((1.-alf)*exp(-argum)+alf*exp(-ho))
        else
          plm(i)=plm(i)*exp(-argum)
        end if
        i=i+1
        go to 1
    2 continue
        hi=tkg/am(4)
        argum=rh/hi*hr
        if(argum.ge.100.)argum=100.
        plm(12)=plm(12)*exp(-argum)
      ! hot O 
      tko=bolc*plm(14)/ge/at1
      hoh=tko/am(3)
      arg=rh/hoh*hr
      if(arg.ge.100.)arg=100.
      plm(13)=plm(13)*exp(-arg)
      return
      end
