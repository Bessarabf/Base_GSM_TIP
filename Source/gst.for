      function gst(alt,tet,vdv,vdu,dolm)
      double precision c,cs,t
      data re/6371.02e5/,g0/980.665/,om/7.272205e-5/
      t=tet
      c=dcos(t)
      cs=c*c
      r=1./(re+alt)
      b=2.*vdv*omt(tet,dolm)
cccc  b=0.
      a=omr(tet,dolm)
      a=(a*a-om*om)/r
      d=1./(1.+3.*cs)
      ds=sqrt(d)
      d=vdv*vdv*.5+(1.+cs)*d*vdu*vdu
cccc  d=0.
c     gst=2.*c*ds*(g0*re*r*re-3.*d)*r
      gst=2.*c*ds*((g0*re*r*re-3.*d)*r+a+b)
      return
      end
