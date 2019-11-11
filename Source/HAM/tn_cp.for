! ver 20.09.19  left and right  first derivative in 1 point is equal
       subroutine tn_ic(an6,vi,rads,n,n1,n2,dt,n0)
c     . . . циклическая прогонка вдоль меридиана
c     . . .    nm=n1+n1-2


       dimension an6(n1,n2,n),vi(n1,n2,n),rads(n)
       allocatable tn(:),v(:),a(:)
     *,            b(:),c(:),f(:)
       data pi/3.1415926/,re/6.371e8/,par/1./
       nm=n1+n1-2
       allocate (tn(nm),v(nm),a(nm)
     *,          b(nm),c(nm),f(nm))

       ns=n1-1
       n_d=n2/2
       dtet=pi/ns
       ot=dt/dtet
       do k=2,n0
        rk=rads(k)+re
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
              a(i)=vp*ot*par/rk
              c(i)=-vm*ot*par/rk
              b(i)=1.+a(i)+c(i)
              f(i)=tn(i)-vm*(1.-par)*ot/rk*(tn(ip)-tn(i))-
     *                   vp*(1.-par)*ot/rk*(tn(i)-tn(im))
            end do
            call cyclp(a,b,c,f,tn,nm)
            do i=1,n1
              an6(i,j,k)=tn(i)
            end do
            do i=n1+1,nm
              an6(i-n1+1,j+n_d,k)=tn(i)
            end do
        end do
       end do
       deallocate (tn,v,a
     *,          b,c,f)
       return
       end

       subroutine tn_jc(an6,vj,rads,n,n1,n2,dt,n0)
c     . . . циклическая прогонка
  
       dimension an6(n1,n2,n),vj(n1,n2,n),rads(n)
       allocatable tn(:),a(:),b(:),c(:),f(:)
       allocate (tn(n2),a(n2),b(n2),c(n2),f(n2)) 
       data pi/3.1415926/,re/6.371e8/,par/1./
       nm=n2
       ns=n1-1
       dfi=2.*pi/n2
       ot=dt/dfi
       do k=2,n0
        rk=rads(k)+re
        do i=2,ns
          tet=pi*(i-1)/ns
          del=rk*sin(tet)
          do j=1,n2
           tn(j)=an6(i,j,k)
          end do
          do j=1,n2
            jm=j-1
            jp=j+1
            if(j.eq.1) jm=n2
            if(j.eq.n2) jp=1
            vm=(vj(i,j,k)-abs(vj(i,j,k)))*.5
            vp=(vj(i,j,k)+abs(vj(i,j,k)))*.5
            a(j)=vp*ot*par/del
            c(j)=-vm*ot*par/del
            b(j)=1.+a(j)+c(j)
            f(j)=tn(j)-vm*(1.-par)*ot/del*(tn(jp)-tn(j))-
     *                 vp*(1.-par)*ot/del*(tn(j)-tn(jm))
          end do
          call cyclp(a,b,c,f,tn,n2)
          do j=2,n2
            an6(i,j,k)=tn(j)
          end do
!!!     left and right  first derivative in 1 point is equal
          an6(i,1,k)=(tn(2)+tn(n2))*.5
        end do
       end do
       deallocate (tn,a,b,c,f)
       return
       end
c
c . . . циклическая прогонка
c . . . a, b, c, f - к-ты трехточки
c . . . a(i)y(i-1)-b(i)y(i)+c(i)y(i+1)=-f(i)
      subroutine cyclp(a,b,c,f,y,nx)
c . . . nx1=nx+1
      dimension a(nx),b(nx),c(nx),f(nx),y(nx)
      allocatable alf(:),bet(:),gam(:)
     *         ,p(:),q(:)
      allocate (alf(nx+1),bet(nx+1),gam(nx+1)
     *         ,p(nx+1),q(nx+1))
       alf(2)=c(2)/b(1)
       gam(2)=a(1)/b(1)
       bet(2)=f(1)/b(1)
       do i=2,nx
         ab=b(i)-a(i)*alf(i)
         alf(i+1)=c(i)/ab
         gam(i+1)=a(i)*gam(i)/ab
         bet(i+1)=(a(i)*bet(i)+f(i))/ab
       end do
       p(nx-1)=bet(nx)
       q(nx-1)=alf(nx)+gam(nx)
       do j=nx-2,1,-1
         p(j)=alf(j+1)*p(j+1)+bet(j+1)
         q(j)=alf(j+1)*q(j+1)+gam(j+1)
       end do
       del=1-alf(nx+1)*q(1)-gam(nx+1)
       y(nx)=(bet(nx+1)+alf(nx+1)*p(1))/del
       do i=1,nx-1
        y(i)=p(i)+y(nx)*q(i)
       end do
       deallocate(alf,bet,gam,p,q)
       return
       end
