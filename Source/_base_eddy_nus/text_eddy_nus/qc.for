      function qc(vic,lam,teta)
      real lam,lamg,lmn,lms
      data rad/57.2957/,vc/3.0e8/
      data fmdn/15./,fmnn/25./,dfn/10./,lmn/45./,dln/25.0/ ! north hemisphere
c     data fmdn/15./,fmnn/25./,dfn/10./,lmn/45./,dln/35.0/ 
      data fmds/15./,fmns/25./,dfs/10./,lms/00./,dls/1.e6/ ! south hemisphere 
c
      lamg=lam*rad
c     !for lmn=0 only!
c     if(lamg.ge.180.) lamg=360.-lamg
      tetag=teta
c
      if(tetag.ge.90.) then
        fms=(fmds+fmns)/2.+cos(lam)*(fmds-fmns)/2.
c       r=((tetag-fms)/dfs)**2
c       r=((tetag-180.+fms)/dfs)**2
        r=((tetag-175.+fms)/dfs)**2
        r=r+((lamg-lms)/dls)**2
c       qc=vic*exp(-r)*2.
c       qc=vic*exp(-r)*0.25
        qc=vic*exp(-r)*0.5
      else
        fmn=(fmdn+fmnn)/2.+cos(lam)*(fmdn-fmnn)/2.
        r=((tetag-fmn)/dfn)**2
        if(lamg.ge.(180.+lmn).and.lmn.le.180.)
     *   lamg=lamg-360.
        r=r+((lamg-lmn)/dln)**2
        qc=vic*exp(-r)
c       if(tetag.le.10.) qc=qc+vc
       end if
c
c     qc=vic*exp(-r)

c     r2=((tetag-10.)/10.)**2
c     if(tetag.ge.90.)r2=((tetag-170.)/10.)**2
c     qc=vic*exp(-r)+vc*exp(-r2)

  900 format(' ',10g12.4)
c     print 900,vic,lam,teta,qc
      return
      end
