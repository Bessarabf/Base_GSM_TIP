      function qf(vic,lam,teta)
      real lam,lamg,lm,lmn,lms
      data rad/57.2957/,vc/0.0e8/
      data fmdn/15./,fmnn/25./,dfn/10./,lmn/195./,dln/60./
      data fmds/0./,fmns/0./,dfs/20./,lms/120./,dls/25./
c
      lamg=lam*rad
c     if(lamg.ge.180.) lamg=360.-lamg !for lm=0 only!
      tetag=teta
c
      if(tetag.ge.30.and.tetag.le.150.) then
        r=0.
        qf=vic*exp(-r)
        goto 1
        else
        qf=0.
        goto 1
      end if
c     if(tetag.ge.90.) then
c       fms=(fmds+fmns)/2.+cos(lam)*(fmds-fmns)/2.
c   !!! precipitation in South Geomagn. Anomaly!!!
c       r=((tetag-125.)/dfs)**2
cc      r=((tetag-175.+fms)/dfs)**2
c       r=r+((lamg-lms)/dls)**2
c       qf=vic*exp(-r)
c     else
c       if(lamg.ge.(180.+lmn).and.lmn.le.180.)
c    *   lamg=lamg-360.
c       if(lamg.le.(lmn-180.).and.lmn.ge.180.)
c    *   lamg=360.+lamg
c       fmn=(fmdn+fmnn)/2.+cos(lam)*(fmdn-fmnn)/2.
c       if(lamg.le.(lmn-180.).and.lmn.ge.180.)
c    *   lamg=360.+lamg
c       r=((tetag-fmn)/dfn)**2
c       r=r+((lamg-lmn)/dln)**2
c       qf=vic*exp(-r)
c      end if
c     qf=vic*exp(-r)+vc*exp(-r2)
    1 continue
  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qc
      return
      end
