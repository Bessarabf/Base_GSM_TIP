      function be(ne,dt,r,rps,vio,vih,vihe,s1,s2)
      dimension c(3)
      data c/5.2e6,8.31e7,2.08e7/
      goto(1,2,3),ne
    1 continue
c
        a=(s1*vih+s2*vihe)*rps
        goto4
    2 continue
c
        a=(s1*vio+s2*vihe)*rps
        goto4
    3 continue
c
        a=(s1*vio+s2*vih)*rps
    4 continue
      be=c(ne)*dt-r-a
      return
      end
