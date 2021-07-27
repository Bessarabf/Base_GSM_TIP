c       dvn :   lambet, alga
c
      function dvn(alt,tet,vdu)
      double precision s,c,cs,d,t
      data re/6371.02e5/
      t=tet
      s=dsin(t)
      c=dcos(t)
      cs=c*c
      d=dsqrt(1.d0+3.d0*cs)
      d=d**3
      dvn=-6.*s*(1.d0+cs)*vdu/((re+alt)*d)
      return
      end
