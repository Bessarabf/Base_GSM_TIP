      subroutine facef(nc,na,alt,dtets,par,kpar,tet1,tet2,tet3,bmpy,
     *           tau,nvg,fac1,fac2,fac3,mast)
      dimension alt(na),par(kpar,na,nc),mast(40)

	allocatable sp(:)
      data e /1.60219e-20/,oe/2.07e7/,oi/3.1e2/,g10/-.30356/,
     *     re/6371.02e5/,ci1/4.23e-10/,ci2/4.28e-10/,ci3/2.58e-10/,
     *     ce1/1.82e-10/,ce11/3.6e-2/,ce2/2.33e-11/,ce21/1.21e-4/,
     *     ce3/2.8e-10/,pi/3.14159265359/,
     *     alf0/4.2e-7/,t0/300./
	 
	allocate (sp(na))

      tet1s=180.-tet1
      tet2s=180.-tet2
      tet3s=180.-tet3
      alf=alf0*t0*2.
      do10i=1,nc
        tets=(i-1)*dtets
        if(tets.eq.tet1.or.tets.eq.tet1s)goto6
          if(tets.eq.tet2.or.tets.eq.tet2s)goto5
            if(tets.eq.tet3.or.tets.eq.tet3s)goto1
              goto10
    1       continue
            if(bmpy.eq.0.)goto10
c             tau1=10./12.*pi
              tau2=pi
c             tau3=14./12.*pi
              tau1=mast(22)/12.*pi
              tau3=mast(30)/12.*pi
              if(bmpy.lt.0.)goto4
                if(tau.lt.tau1.or.tau.gt.tau3)goto10
                  if(tau.eq.tau2)goto10
                    if(tau.gt.tau2)goto3
    2                 continue
                      fac=fac3
                      if(tets.eq.tet3)fac=-fac
                      goto7
    3               continue
                    fac=fac3
                    if(tets.eq.tet3s)fac=-fac
                    goto7
    4             continue
                  if(tau.lt.tau1.or.tau.gt.tau3)goto10
                    if(tau.eq.tau2)goto10
                      if(tau.gt.tau2)goto2
                        goto3
    5     continue
          fac=-fac2*sin(tau)
          goto7
    6   continue
        fac=fac1*sin(tau)
        if(nvg.eq.0)goto7
          fac=0.
    7   continue
        tet=tets/180.*pi
        ct=cos(tet)
        sk=sqrt(1.+3.*ct*ct)
        h=alt(1)
        sip=0.
        do8j=1,na
          roq=(re/(re+alt(j)))**3
          b=-g10*roq*sk
          ee=e/b
          ome=oe*b
          omi=oi*b
          cm=par(4,j,i)
          co2=par(1,j,i)
          cn2=par(2,j,i)
          co=par(3,j,i)
c         tn=par(5,j,i)
          te=par(5,j,i)
c         b=sqrt(tn)
          b=sqrt(te)
          fi=ci1*co2+ci2*cn2+ci3*co
c         fe=ce1*(1.+ce11* b)*b*co2+ce2*(1.-ce21*tn)*tn*cn2+ce3*b*co
          fe=ce1*(1.+ce11*b)*b*co2+ce2*(1.-ce21*te)*te*cn2+ce3*b*co
          bu=ome*ome
          cu=omi*omi
          co2=fi*fi
          cn2=fe*fe
          b=1./(bu+cn2)
          co=1./(cu+co2)
          sp(j)=cm*ee*(omi*fi*co+ome*fe*b)
          if(j.eq.1)goto8
            hp=alt(j)
            sip=(sp(j)+sp(j-1))*.5*(hp-h)
            h=hp
    8   continue
        do9j=1,na
          cm=par(4,j,i)
c         tn=par(5,j,i)
          te=par(5,j,i)
c         b=cm-sp(j)/(e*sip)*fac*tn/(alf*cm)
          b=cm-sp(j)/(e*sip)*fac*te/(alf*cm)
          if(b.lt.1.)b=1.
          par(4,j,i)=b
    9   continue
   10 continue
      deallocate (sp)
      return
      end
