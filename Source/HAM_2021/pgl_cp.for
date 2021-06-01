       subroutine pgl_ic(pgl,rads,kpars,nh,its,ids,dt,n0,np)
c     . . . циклическая прогонка вдоль меридиана
c     . . .    nm=its+its-2


       dimension pgl(kpars,nh,its,ids),rads(nh)
       allocatable tn(:),v(:),a(:)
     *,            b(:),c(:),f(:)
       data pi/3.1415926/,re/6.371e8/,par/1./
      
       nm=its+its-2
       allocate( tn(nm),v(nm),a(nm)
     *,            b(nm),c(nm),f(nm))
       ns=its-1
       n_d=ids/2
       dtet=pi/ns
       ot=dt/dtet
       do k=2,n0
        rk=rads(k)+re
        do j=1,n_d
c . . . переход к меридиональному кругу
            do i=1,its
              v(i)=pgl(11,k,i,j)
              tn(i)=pgl(np,k,i,j)
            end do
            do i=its+1,nm
              v(i)=-pgl(11,k,i-its+1,j+n_d)
              tn(i)=pgl(np,k,i-its+1,j+n_d)
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
            do i=1,its
              pgl(np,k,i,j)=tn(i)
            end do
            do i=its+1,nm
              pgl(np,k,i-its+1,j+n_d)=tn(i)
            end do
        end do
       end do
       deallocate( tn,v,a,b,c,f)
       return
       end
      
	subroutine pgl_jc(pgl,rads,kpars,nh,its,ids,dt,n0,np)
c     . . . циклическая прогонка

       dimension pgl(kpars,nh,its,ids),rads(nh)
       allocatable tn(:),a(:),b(:),c(:),f(:)
       data pi/3.1415926/,re/6.371e8/,par/1./
       allocate (tn(ids),a(ids),b(ids),c(ids),f(ids))
       nm=ids
       ns=its-1
       dfi=2.*pi/ids
       ot=dt/dfi
       do k=2,n0
        rk=rads(k)+re
        do i=2,ns
          tet=pi*(i-1)/ns
          del=rk*sin(tet)
          do j=1,ids
           tn(j)=pgl(np,k,i,j)
          end do
          do j=1,ids
            jm=j-1
            jp=j+1
            if(j.eq.1) jm=ids
            if(j.eq.ids) jp=1
            vm=(pgl(12,k,i,j)-abs(pgl(12,k,i,j)))*.5
            vp=(pgl(12,k,i,j)+abs(pgl(12,k,i,j)))*.5
            a(j)=vp*ot*par/del
            c(j)=-vm*ot*par/del
            b(j)=1.+a(j)+c(j)
            f(j)=tn(j)-vm*(1.-par)*ot/del*(tn(jp)-tn(j))-
     *                 vp*(1.-par)*ot/del*(tn(j)-tn(jm))
          end do
          call cyclp(a,b,c,f,tn,ids)
          do j=2,ids
            pgl(np,k,i,j)=tn(j)
          end do
        end do
!!!     left and right  first derivative in 1 point is equal
        pgl(np,k,i,1)=(tn(2)+tn(ids))*.5

       end do
       deallocate(tn,a,b,c,f)
       return
       end
