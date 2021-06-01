      function qp(vic,lam,teta)
      real lam,lamg,lm,lmn,lms
      data rad/57.2957/,vc/0.0e8/
      data fmdn/15./,fmnn/20./,dfn/10./,lmn/195./,dln/60./
c     data fmdn/15./,fmnn/25./,dfn/10./,lmn/195./,dln/60./
c     data fmdn/15./,fmnn/25./,dfn/10./,lmn/165./,dln/30./
c     data fmdn/15./,fmnn/30./,dfn/10./,lmn/165./,dln/60./
c     data fmds/0./,fmns/0./,dfs/25./,lms/150./,dls/1.e6/
c     data fmds/0./,fmns/0./,dfs/20./,lms/150./,dls/25./
      data fmds/0./,fmns/0./,dfs/20./,lms/120./,dls/25./
c
      lamg=lam*rad
c     if(lamg.ge.180.) lamg=360.-lamg !for lm=0 only!
      tetag=teta
c
      if(tetag.ge.90.) then
        fms=(fmds+fmns)/2.+cos(lam)*(fmds-fmns)/2.
c   !!! precipitation in South Geomagn. Anomaly!!!
        r=((tetag-125.)/dfs)**2
cc      r=((tetag-175.+fms)/dfs)**2
        r=r+((lamg-lms)/dls)**2
        qp=vic*exp(-r)
      else
c       if(lamg.ge.(180.+lmn).and.lmn.le.180.)
c    *   lamg=lamg-360.
c       if(lamg.le.(lmn-180.).and.lmn.ge.180.)
c    *   lamg=360.+lamg
        fmn=(fmdn+fmnn)/2.+cos(lam)*(fmdn-fmnn)/2.
        if(lamg.le.(lmn-180.).and.lmn.ge.180.)
     *   lamg=360.+lamg
        r=((tetag-fmn)/dfn)**2
        r=r+((lamg-lmn)/dln)**2
c       qp=vic*exp(-r)*0.40
c       qp=vic*exp(-r)*0.30
        qp=vic*exp(-r)*0.15
       end if
c     qp=vic*exp(-r)+vc*exp(-r2)

  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qc
      return
      end
