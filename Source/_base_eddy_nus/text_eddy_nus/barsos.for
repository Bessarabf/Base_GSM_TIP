      subroutine barsos(an,an6,rp,g,am,n,n1,n2,l)
      dimension an(n1,n2,n)
     *         ,an6(n1,n2,n),rp(n),g(n)
      data bk/1.38e-16/,ves/0.5/
      f2=bk/am
      do 1 i=1,n1
       do 2 j=1,n2
        do 3 k=l,n
          tn=an6(i,j,k-1)
          tv=an6(i,j,k)
          h1=f2*tn/g(k-1)
          h2=f2*tv/g(k)
          alf=(ves/h1+(1.-ves)/h2)*rp(k-1)
          ss=alog(an(i,j,k-1)*tn/tv)-alf
          an(i,j,k)=exp(ss)
    3   continue
    2  continue
    1 continue
      return
      end
