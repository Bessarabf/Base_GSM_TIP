      real function nuOhot(ne,ti,tOhot)
      data s1/.37/,s2/.38/
        t=ti+tOhot
        goto(2,3),ne
    2   continue
c
          a=t**s1
          goto4
    3   continue
c
          a=t**s2
    4   continue
        nuOhot=a
        return
      end
