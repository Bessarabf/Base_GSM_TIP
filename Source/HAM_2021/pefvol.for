      function pefvol(delta,utt,pkp,bmpy,alt,tet,fi,j,j1,nl2)
      real l,l0,l1,l3
      double precision pi
c     data fic1/8.6e11/,fic2/3.2e12/,fip/9.e12/,pi/3.14159265359d0/,
      data fic1/8.6e11/,fic2/4.0e12/,fip/9.e12/,pi/3.14159265359d0/,
     *l0/14.9/,l1/8.5/,l3/33.2/,re/6371.02e5/,
     *c1/1.908395e-1/,c2/4.6479083e-1/,c3/9.184622e-1/,
     *c4/2.029585/,c5/4.4961194e-2/
      pp=pi+pi
c     print 113
  113 format(' pefvol')
      if(j.eq.1.or.j.eq.nl2)go to 3
        l=(re+alt)/(re*sin(tet)**2)
c       tau = fi-fih(delta,utt)+pi
        call magsm(utt,delta,fi,phism,1)
        tau=pi+phism
        if(tau.ge.pp)tau=tau-pp
        if(tau.lt.0.)tau=pp+tau
        if(l.gt.l0)go to 3
          if(pkp.gt.2..and.pkp.le.4.)goto1
            pefvol=fic1*(l/l0)**2*sin(tau)
            return
    1     continue
          if(l.gt.l1)goto2
            tau=tau-c1
            pefvol=fic2*c2*(l/l1)**2*sin(tau)
            return
    2     continue
            tau=tau-c3+.34*alog(l)
            pefvol=fic2*(l/l0)**1.365*sin(tau)
            return
    3 continue
        if(bmpy.ne.0.)goto4
          pefvol=0.
          if(j.eq.1.or.j.eq.nl2)return
          c=0.
          go to 5
    4   continue
          c=fip
          if(bmpy.gt.0..and.j.gt.j1.or.bmpy.lt.0..and.j.lt.j1)c=-fip
          if(j.ne.1.and.j.ne.nl2)go to 5
            pefvol=c*c5
            return
    5     continue
          d=sqrt(l0/l)*sin(tau)
          if(pkp.gt.2..and.pkp.le.4.)goto6
            d=fic1*d
            goto7
    6     continue
            d=fic2*d
    7     continue
          if(l.gt.l3)goto8
            pefvol=c*c4*(sqrt(l/l0)-1.)**2/l+d
            return
    8     continue
            pefvol=c*(c5-1./l)+d
            return
      end

