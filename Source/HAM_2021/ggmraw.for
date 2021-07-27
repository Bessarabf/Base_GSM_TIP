      subroutine ggmraw(art,dolg,tet,dolm,tetm)
! art - var of convertation 
! 0 - geo to mag coordinate
! 1 - mag to geo coordinate 
! input and output in gratitudes
      integer art
      double precision
     *    zpi,faktor,cbg,ci,si,xlm,bm,cbm,sbm,
     *    clm,slm,sbg,bg,slg,clg,xlg,ylg
      zpi=6.28318530718
      faktor=0.01745329252
      cbg=11.4*faktor
      ci=dcos(cbg)
      si=dsin(cbg)
      if(art.eq.0) go to 10
        xlm=dolm
        bm=90.-tetm
        cbm=dcos(bm*faktor)
        sbm=dsin(bm*faktor)
        clm=dcos(xlm*faktor)
        slm=dsin(xlm*faktor)
        sbg=sbm*ci-cbm*clm*si
        bg=dasin(sbg)
        cbg=dcos(bg)
        slg=(cbm*slm)/cbg
        clg=(sbm*si+cbm*clm*ci)/cbg
        if(clg.gt.1..or.(1.-clg).lt.1.e-10)clg=1.
        if(clg.lt.-1..or.(1.+clg).lt.1.e-10)clg=-1.
        xlg=dacos(clg)
        if(slg.lt.0.0) xlg=zpi-dacos(clg)
        bg=bg/faktor
        xlg=xlg/faktor
        xlg=xlg-69.8
        if(xlg.lt.0.0) xlg=xlg+360.0
        tet=90.-bg
        dolg=xlg
        go to 20
   10   bg=90.-tet
        xlg=dolg
        ylg=xlg+69.8
        cbg=dcos(bg*faktor)
        sbg=dsin(bg*faktor)
        clg=dcos(ylg*faktor)
        slg=dsin(ylg*faktor)
        sbm=sbg*ci+cbg*clg*si
        bm=dasin(sbm)
        cbm=dcos(bm)
        slm=(cbg*slg)/cbm
        clm=(-sbg*si+cbg*clg*ci)/cbm
        if(clm.lt.-1.0) clm=-1.0 
	  if(clm.gt.1.0) clm=1.0 
	  xlm=dacos(clm)
        if(slm.lt.0.0) xlm=zpi-dacos(clm)
        bm=bm/faktor
        xlm=xlm/faktor
        dolm=xlm
        if(abs(tet-180.).lt.1.e-3)dolm=0.
        tetm=90.-bm
   20 continue
      return
      end
