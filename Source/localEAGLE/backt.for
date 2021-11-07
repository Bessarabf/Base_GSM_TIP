      subroutine backt(i1,i2,delu,bet,gam)

      dimension delu(*),bet(*),gam(*)
      i4=i2-1
      do1i=i1,i4
        m=i4-i+i1
        mp=m+1
        delu(m)=bet(m)*delu(mp)+gam(m)
    1 continue
      return
      end
