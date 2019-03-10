      subroutine g52041(ni,nim,nj,x,y,jww,psdv)
      dimension x(nim),y(nim),psdv(nim,nim),pr(50,50)
      jww=52041
      mi=ni+1
      mj=nj+1
      do 2 j=1,nj
        do 1 i=1,ni
          pr(i,j)=psdv(i,j)
    1   continue
        pr(mi,j)=pr(1,j)
    2 continue
      do 3 i=1,ni
        pr(i,mj)=pr(i,1)
    3 continue
      pr(mi,mj)=pr(2,2)
      do 5 j=2,nj-1
        tg=y(j)
        do 4 i=1,ni
          dg=x(i)
          psdv(i,j)=g52051(dg,ni,nim,nj,pr,tg,x,y)
    4   continue
    5 continue
      tg=0.
      dg=30.
      p=g52051(dg,ni,nim,nj,pr,tg,x,y)
      do 6 i=1,ni
        psdv(i,1)=p
    6 continue
      tg=180.
      p=g52051(dg,ni,nim,nj,pr,tg,x,y)
      do 7 i=1,ni
        psdv(i,nj)=p
    7 continue
      return
      end
      function g52051(dm,ni,nim,nj,pr,tm,x,y)
      dimension pr(nim,nim),x(nim),y(nim)
      ks=1
      call g52t61(ks,dg,tg,dm,tm)
      i=j52t62(nim,ni,dg,x)
      j=j52t62(nim,nj,tg,y)
      k=i+1
      l=j+1
      r=x(i)
      cx=(dg-r)/(x(k)-r)
      r=y(j)
      cy=(tg-r)/(y(l)-r)
      r=1.-cy
      s=1.-cx
      a=r*s
      b=r*cx
      c=cy*s
      d=cy*cx
      r=pr(i,j)*a+pr(k,j)*b
      s=pr(i,l)*c+pr(k,l)*d
      g52051=r+s
      return
      end

