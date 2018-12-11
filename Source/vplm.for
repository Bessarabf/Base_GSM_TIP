      subroutine vplm(p1,p2,tet)
      double precision f,tr,sn,cs,tn2,si
      f=0.01745329252
      tr=(90.-tet)*f
      sn=dsin(tr)
      cs=dcos(tr)
      tn2=(sn/cs)*2
      si=datan(tn2)
      sn=dsin(si)
      cs=dcos(si)
      a1=p1
      a2=p2
      p1=-(a1*sn+a2*cs)
      p2=a2*sn-a1*cs
      return
      end
