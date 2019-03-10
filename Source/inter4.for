      subroutine inter4(x1,x2,xi,vn1,vn2,plm,kpart,log)
      dimension plm(kpart),vn1(kpart),vn2(kpart)

      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
c
      if(log.eq.0)go to 11
        do 1 i=1,3
          p1=alog(vn1(i))
          p2=alog(vn2(i))
          plm(i)=exp(w1*p2+w2*p1)
    1   continue
        do 2 i=4,8
          plm(i)=w1*vn2(i)+w2*vn1(i)
    2   continue
        go to 12
   11 continue
        do 13 i=1,8
          plm(i)=w1*vn2(i)+w2*vn1(i)
   13   continue
   12 continue
    
      return
      end
