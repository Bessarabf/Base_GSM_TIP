      subroutine inter2(x1,x2,xi,vn1,vn2,plm)
      dimension plm(6),vn1(6),vn2(6)
      a1=amin1(x1,x2)
      a2=amax1(x1,x2)
      ! 10.03.2018 first and lust tube is vertical
      if(xi.lt.a1.and.a1.eq.x1.or.xi.gt.a2.and.a2.eq.x1)then
        do i=1,6
          plm(i)=vn1(i)
        end do
      end if
      if(xi.lt.a1.and.a1.eq.x2.or.xi.gt.a2.and.a2.eq.x2)then
        do i=1,6
          plm(i)=vn2(i)
        end do
      end if
      if(xi.ge.a1.and.xi.le.a2)then
        hi=x2-x1
        w1=(xi-x1)/hi
        if(xi.lt.a1.or.xi.gt.a2)print*,' ints x1,xi,x2  ',x1,xi,x2
        w2=1.-w1
c       if(vn1(1).le.0.)vn1(1)=1.e15
c       if(vn2(1).le.0.)vn2(1)=1.e15
c       p1=alog(vn1(1))
c       p2=alog(vn2(1))
c       plm(1)=exp(w1*p2+w2*p1)
        do i=1,6
          plm(i)=w1*vn2(i)+w2*vn1(i)
        end do
      end if
      return
      end
