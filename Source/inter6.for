      subroutine inter6(x1,x2,xi,vn1,vn2,kpart,log,nv,l,pm)
      dimension pm(kpart,nv),vn1(kpart),vn2(kpart)
      data al/0.54e-10/
      hi=x2-x1
      w1=(xi-x1)/hi
      w2=1-w1
      do 1 i=1,3
        if(vn1(i).lt.al)vn1(i)=al
        if(vn2(i).lt.al)vn2(i)=al
        p1=vn1(i)
        p2=vn2(i)
        if(log.eq.0)go to 7
          p1=alog(vn1(i))
          p2=alog(vn2(i))
    7   continue
c       if(l.eq.5)print 77,l,nv,log,i,p1,p2
  77  format (' ',6g12.4)
        pm(i,l)=w1*p2+w2*p1
    1 continue
      do 2 i=4,8
        pm(i,l)=w1*vn2(i)+w2*vn1(i)
    2 continue
      return
      end
