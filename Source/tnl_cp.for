       subroutine tnl_ic(an6,an1,an2,an3,vi,ctd,rads,n,n1,n2,dt,n0)
c     . . . теплопроводность вдоль меридиана и долготы
c     . . . циклическая прогонка вдоль меридиана
c     . . .    nm=n1+n1-2
       
      
       dimension an6(n1,n2,n),an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n)
     *,          vi(n1,n2,n),ctd(n),rads(n)
     
       allocatable tn(:),v(:),a(:),b(:),c(:),f(:),alyam(:)
       data pi/3.1415926/,re/6.371e8/,bk/1.38e-16/,par/1./
       nm=n1+n1-2
       allocate (tn(nm),v(nm),a(nm),b(nm),c(nm),f(nm)
     *,          alyam(nm))
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
c . . . теплопроводность
              al=18.6*tn(i)**0.84*an1(i,j,k)
              al=al+27.2*tn(i)**0.8*an2(i,j,k)
              al=al+67.1*an6(i,j,k)**0.71*an3(i,j,k)
              alyam(i)=al/(an1(i,j,k)+an2(i,j,k)+
     *                  an3(i,j,k))
c         . . . + турбулентная теплопроводность
              rocp=(3.5*(an1(i,j,k)+an2(i,j,k))+2.5*
     *              an3(i,j,k))*bk
              alyam(i)=alyam(i)+rocp*ctd(k)
            end do
            do i=n1+1,nm
              ipol=i-n1+1
              v(i)=-vi(ipol,j+n_d,k)
              tn(i)=an6(ipol,j+n_d,k)
              al=18.6*tn(i)**0.84*an1(ipol,j,k)
              al=al+27.2*tn(i)**0.8*an2(ipol,j,k)
              al=al+67.1*an6(ipol,j,k)**0.71*an3(ipol,j,k)
              alyam(i)=al/(an1(ipol,j,k)+an2(ipol,j,k)+
     *                  an3(ipol,j,k))
              rocp=(3.5*(an1(ipol,j,k)+an2(ipol,j,k))+2.5*
     *              an3(ipol,j,k))*bk
              alyam(i)=alyam(i)+rocp*ctd(k)
            end do
            do i=1,nm
              im=i-1
              ip=i+1
              if(i.eq.nm) ip=1
              if(i.eq.1) im=nm
c . . . новая переменная по широте
              tetn=dtet*(i-1)
              sin_t=abs(sin(tetn))
              if(sin_t.lt.0.03) sin_t=0.03
              del=rk*sin_t*dtet
              vm=(v(i)-abs(v(i)))*.5
              vp=(v(i)+abs(v(i)))*.5
              arg=tetn+dtet*.5
              ap_2=abs(sin(arg))*(alyam(ip)+alyam(i))*.5
              arg=tetn-dtet*.5
              am_2=abs(sin(arg))*(alyam(i)+alyam(im))*.5
              a(i)=(vp*par+am_2/del)*ot/rk
              c(i)=(-vm*par+ap_2/del)*ot/rk
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
       deallocate (tn,v,a,b,c,f,alyam)
       return
       end


       subroutine tnl_jc(an6,an1,an2,an3,vj,ctd,rads,n,n1,n2,dt,n0)
c     . . . циклическая прогонка
       
       
       dimension an6(n1,n2,n),an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n)
     *,          vj(n1,n2,n),ctd(n),rads(n)
     
       allocatable tn(:),a(:),b(:),c(:),f(:),alyam(:)
       allocate (tn(n2),a(n2),b(n2),c(n2),f(n2),alyam(n2))
       data pi/3.14159/,re/6.371e8/,bk/1.38e-16/,par/1./
       ns=n1-1
       dfi=2.*pi/n2
       dtet=pi/ns
       ot=dt/dfi
       do k=2,n0
        rk=rads(k)+re
        do i=2,ns
          tet=dtet*(i-1)
          del=rk*sin(tet)
          do j=1,n2
           tn(j)=an6(i,j,k)
c . . . теплопроводность
           al=18.6*tn(j)**0.84*an1(i,j,k)
           al=al+27.2*tn(j)**0.8*an2(i,j,k)
           al=al+67.1*an6(i,j,k)**0.71*an3(i,j,k)
           alyam(j)=al/(an1(i,j,k)+an2(i,j,k)+
     *              an3(i,j,k))
c         . . . + турбулентная теплопроводность
              rocp=(3.5*(an1(i,j,k)+an2(i,j,k))+2.5*
     *              an3(i,j,k))*bk
              alyam(j)=alyam(j)+rocp*ctd(k)
          end do
          do j=1,n2
            jm=j-1
            jp=j+1
            if(j.eq.1) jm=n2
            if(j.eq.n2) jp=1
            vm=(vj(i,j,k)-abs(vj(i,j,k)))*.5
            vp=(vj(i,j,k)+abs(vj(i,j,k)))*.5
            ap_2=(alyam(jp)+alyam(j))*.5
            am_2=(alyam(j)+alyam(jm))*.5
            a(j)=(vp*par+am_2/(del*dfi))/del*ot
            c(j)=(-vm*par+ap_2/(del*dfi))/del*ot
            b(j)=1.+a(j)+c(j)
            f(j)=tn(j)-vm*(1.-par)*ot/del*(tn(jp)-tn(j))-
     *                 vp*(1.-par)*ot/del*(tn(j)-tn(jm))
          end do
          call cyclp(a,b,c,f,tn,n2)
          do j=1,n2
            an6(i,j,k)=tn(j)
          end do
        end do
       end do
       deallocate (tn,a,b,c,f,alyam)
       return
       end

