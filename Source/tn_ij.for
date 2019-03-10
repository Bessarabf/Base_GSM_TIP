C
      subroutine tngoim(an6,vi,rads,n,n1,n2,dt,n0)
     
      dimension an6(n1,n2,n),vi(n1,n2,n),rads(n)
      allocatable tn(:)
      allocate (tn(n1))
       data pi/3.1415926/,re/6.371e8/
       ns=n1-1
       nm=n2-1
       df=pi/ns
       do1k=1,n0
        rk=rads(k)+re
        do 2 j=1,n2
         do 3 i=2,ns
          am=(vi(i,j,k)-abs(vi(i,j,k)))/(df*rk)*.5
          ap=(vi(i,j,k)+abs(vi(i,j,k)))/(df*rk)*.5
          a=am*(an6(i+1,j,k)-an6(i,j,k))
          b=ap*(an6(i,j,k)-an6(i-1,j,k))
          tn(i)=an6(i,j,k)-dt*(a+b)
  3      continue
         do 4 i=2,ns
          an6(i,j,k)=tn(i)
  4      continue
  2     continue
  1    continue
       deallocate (tn)
       return
       end
C
      subroutine tngojm(an6,vj,rads,n,n1,n2,dt,n0)
      
       dimension an6(n1,n2,n),vj(n1,n2,n),rads(n)
      allocatable tn(:)
      allocate (tn(n2))
       data pi/3.14159/,re/6.371e8/
       ns=n1-1
       da=2.*pi/n2
       do 1 k=2,n0
        rk=rads(k)+re
        do 2 i=2,ns
          fi=pi*(i-1)/ns
          do 3 j=1,n2
           jm=j-1
           jp=j+1
           if(j.eq.1) jm=n2
           if(j.eq.n2) jp=1
           am=(vj(i,j,k)-abs(vj(i,j,k)))/(rk*da*sin(fi))
           ap=(vj(i,j,k)+abs(vj(i,j,k)))/(rk*da*sin(fi))
           a=am*(an6(i,jp,k)-an6(i,j, k))
           b=ap*(an6(i,j,k)-an6(i,jm,k))
           tn(j)=an6(i,j,k)-dt*(a+b)*.5
  3       continue
          do 4 j=1,n2
           an6(i,j,k)=tn(j)
  4       continue
  2      continue
  1    continue
       deallocate (tn)
       return
       end
C
c . . . Перенос через полюс вдоль меридиана
c . . . явная схема
      subroutine tngoim_a(an6,vi,rads,n,n1,n2,dt,n0)
  
      dimension an6(n1,n2,n),vi(n1,n2,n),rads(n)
      allocatable tn(:),tn1(:),v(:)
      
      data pi/3.1415926/,re/6.371e8/
      nm=n1+n1-2
      allocate (tn(nm),tn1(nm),v(nm))
      ns=n1-1
      n_d=n2/2
      dtet=pi/ns
      ot=dt/dtet
      do k=2,n0
       rk=rads(k)+re
       ot=ot/rk
       do j=1,n_d
c . . . переход к меридиональному кругу
           do i=1,n1
             v(i)=vi(i,j,k)
             tn(i)=an6(i,j,k)
           end do
           do i=n1+1,nm
             v(i)=-vi(i-n1+1,j+n_d,k)
             tn(i)=an6(i-n1+1,j+n_d,k)
           end do
           do i=1,nm
             im=i-1
             ip=i+1
             if(i.eq.nm) ip=1
             if(i.eq.1) im=nm
             vm=(v(i)-abs(v(i)))*.5
             vp=(v(i)+abs(v(i)))*.5
             a=vm*(tn(ip)-tn(i))
             b=vp*(tn(i)-tn(im))
             tn1(i)=tn(i)-ot*(a+b)
           end do
c . . . обратный переход
           do i=1,n1
             an6(i,j,k)=tn1(i)
           end do
           do i=n1+1,nm
             an6(i-n1+1,j+n_d,k)=tn1(i)
           end do
        end do
       end do
       deallocate (tn,tn1,v)
       return
       end