c   nu:  lambet-tube-trubka-cyclt2
c   nu:  rik-lambet-tube-trubka-cyclt2
c   nu:  pitn-alga-tube-trubka-cyclt2
c   nu:  pqji-piqj- alga-tube-trubka-cyclt2
c   nu:  piqj- alga-tube-trubka-cyclt2
c
      real function nu(ne,ti,tn)
      data s1/.37/,s2/.38/,s3/-1.5/
      goto(1,1,4),ne
    1 continue
        t=ti+tn
        goto(2,3),ne
    2   continue
c
          a=t**s1
          goto5
    3   continue
c
          a=t**s2
          goto5
    4   continue
c		
  	 
          a=ti**s3
    5   continue
        nu=a
        return
      end
