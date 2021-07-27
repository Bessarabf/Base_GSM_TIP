      subroutine gu(ll,rads,nh,par,kpars,ncs,nl,nvar,ae,be,ge)
      dimension rads(nh),par(kpars,nh,ncs),ae(nl),be(nl),ge(nl)
      data re/6371.02e5/,pi/3.14159265359/
      h=rads(1)
      t=pi*.5
      b=bdip(h,t)
      r=re+h
      roq=(re/r)**3
      j=(ncs+1)/2
      cm=par(6,1,j)
      co2=par(1,1,j)
      cn2=par(2,1,j)
      co=par(3,1,j)
      tn=par(7,1,j)
      vnu=-par(10,1,j)
      vnv=par(12,1,j)
      call conducn(b,cm,co2,cn2,co,tn,sp,sh)
c     ae(ll)=sp
      ae(ll)=sh/r
c     be(ll)=sh
      be(ll)=-sp*re/(r*r)
      r=0.
      if(nvar.ne.2)r=b*(sh*vnu+sp*vnv)
      ge(ll)=r
      return
      end
