      function qa(vic,lam,teta)
      real lam,lamg,lm
      data rad/57.2957/,vc/0.0e8/
c     data fmd/15./,fmn/25./,df/5./,lm/210./,dl/40./
      data fmd/15./,fmn/25./,df/5./,lm/210./,dl/60./
c     data fmd/20./,fmn/20./,df/5./,lm/180./,dl/60./
c     data fmd/20./,fmn/20./,df/5./,lm/210./,dl/75./
c
      lamg=lam*rad
      tetag=teta
c
      fm=(fmd+fmn)/2.+cos(lam)*(fmd-fmn)/2.
      r=((tetag-fm)/df)**2
c      if(tetag.ge.90.)  r=((tetag-175.+fm)/df)**2
       if(tetag.ge.90.)  r=((tetag-180.+fm)/df)**2
c       if(lamg.ge.(180.+lm).and.lm.le.180.)
c    *   lamg=lamg-360.
        if(lamg.le.(lm-180.).and.lm.ge.180.)
     *   lamg=360.+lamg
      r=r+((lamg-lm)/dl)**2
      qa=vic*exp(-r)

c     r2=((tetag-10.)/10.)**2
c     if(tetag.ge.90.)r2=((tetag-170.)/10.)**2
c     qa=vic*exp(-r)+vc*exp(-r2)

  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qa
      return
      end
