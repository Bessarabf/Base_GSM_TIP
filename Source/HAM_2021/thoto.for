      subroutine thotO(pgl,pgi,cHot,tHot,
     *                 qHot,rads,kpars,ins,nh,its,ids,dts,mass)
  
      dimension pgl(kpars,nh,its,ids),pgi(ins,nh,its,ids),
     *          tHot(its,ids,nh),cHot(its,ids,nh),
     *          qHot(its,ids,nh),
     *          rads(nh),mass(*)
    
     *,         alyam(3)
      allocatable pa(:),pb(:)
      allocate (pa(NH+5),pb(NH+5))
      data re/6.371e8/,pi/3.1415926/,bk/1.38e-16/
    	!! 
!	 open (444,file='qhot.tmp')
!	rewind(444)
      !!!
	cr=pi/180.
      klik=0
      n0=nh
      n0m=n0-1
      itsm1=its-1
      itsm2=its-2
      dteta=pi/itsm1
      ddol=2.*pi/ids
c
c
      do 1 i=2,itsm1
       teta=dteta*(i-1)
       sin t=sin(teta)
       sin p=sin(teta+dteta)
       sin m=sin(teta-dteta)
       do 2 j=1,ids
c
        pa(2)=0.
        pb(2)=thot(i,j,1)
        jm=j-1
        jp=j+1
        if(j.eq.1) jm=ids
        if(j.eq.ids) jp=1
       
        do 3 k=2,n0m
c
      ! teploprovodnost v 3 tochkah  
        do 20 k1=1,3
          kv=k-2+k1
          alyam(k1)=75.9*thot(i,j,kv)**0.69
   20   continue
c
        
        rc=1.5*cHot(i,j,k)*bk
        abs vr= abs(pgl(10,k,i,j))
        con1=pgl(10,k,i,j)-abs vr
        con2=pgl(10,k,i,j)+abs vr
c
        drad=rads(k+1)-rads(k)
        drad1=rads(k)-rads(k-1)
c
        rk=rads(k)+re
c
        a=((alyam(3)+alyam(2))/drad-rc*con1)/(rads(k+1)-
     *     rads(k-1))
        c=((alyam(2)+alyam(1))/drad1+rc*con2)/(rads(k+1)-
     *     rads(k-1))
        b=a+c+rc/dts
c
        davl0=(pgl(10,k+1,i,j)-pgl(10,k-1,i,j))/(drad+drad1)
!        davl1=(pgl(11,k,i+1,j)*sin p-
!     *         pgl(11,k,i-1,j)*sin m)/dteta*.5
        davl1=(pgl(11,k,i+1,j)-
     *         pgl(11,k,i-1,j))/dteta*.5+pgl(11,k,i,j)*cos(teta)/sin t
        davl2=(pgl(12,k,i,jp)-pgl(12,k,i,jm))/ddol*.5
c        davl=davl0+(davl1+davl2)/(rk*sin t)
        davl=davl0+davl1/rk+davl2/(rk*sin t)

        ha=rc*tHot(i,j,k)/dts
c
!. . . lost neutral      
       !c razn=tHot(i,j,k)-pgl(7,k,i,j)
	 razn=-pgl(7,k,i,j)

    ! proizv=7.8e-9*bk*cHot(i,j,k)*
    ! * 	      (0.061*pgl(3,k,i,j)+0.075*pgl(2,k,i,j))
!!!	 proizv=1.715e-12*bk*razn*cHot(i,j,k)*
!!     *        pgl(3,k,i,j)*sqrt(thot(i,j,k)+pgl(7,k,i,j))	!hard sphere 

	  proizv=1.143e-12*bk*cHot(i,j,k)*
     *    (pgl(3,k,i,j)*sqrt(pgl(7,k,i,j)+pgl(7,k,i,j))+!hard sphere 
     *     pgl(2,k,i,j)*2.031*sqrt(pgl(7,k,i,j)+pgl(7,k,i,j)*0.571))!approx:Brun													!approx: Brunelly, Namgaladze  
         
          qlost=proizv*razn
!. . .  lost for O+
          if(pgi(6,k,i,j).gt.1) then ! эх Ёртэю 0 ёю 175 ъь
           prion=1.676e-10*bk*cHot(i,j,k)*
     *		  (tHot(i,j,k)+pgi(6,k,i,j))**0.37*pgi(1,k,i,j)
           prion=prion*.5 ! ЁрчюсЁрЄ№ё  ё ўрёЄюЄющ! єьхэ№°хэ т 2 Ёрчр 
          else 
            prion=0.
          end if	 
	  qlost=qlost+prion*(-pgi(6,k,i,j))
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !c if(abs(razn).Le.0.1) then
       !c   qlost=0.
       !c else
       !c    qlost=7.8e-9*bk*razn*cHot(i,j,k)*
    !c * 	   (0.061*pgl(3,k,i,j)+0.075*pgl(2,k,i,j))
       !	 qlost=1.715e-12*bk*razn*cHot(i,j,k)*
       !	       pgl(3,k,i,j)*sqrt(thot(i,j,k)+pgl(7,k,i,j))	!hard sphere 
!															!approx: Brunelly, Namgaladze   
!	     qlost=0.
!	if(qhot(i,j,k).lt.qlost) print*,i,j,k,qhot(i,j,k),qlost
      !c  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fk=qHot(i,j,k) + ha-qlost
        b=b+proizv+prion
! . . .  ╙ўхЄ чэрър фртыхэш 
        if(davl.lt.0.) then
          fk=fk-davl*cHot(i,j,k)*bk*Thot(i,j,k)
        else
          b=b+davl*cHot(i,j,k)*bk
        end if
! . . .  ╧Ё ьр  яЁюуюэър
        pa(k+1)=a/(b-pa(k)*c)
        pb(k+1)=(c*pb(k)+fk)/(b-pa(k)*c)
!	  if(razn.lt.0) then
!         
!	    write (444,'(3I4,1p(7E9.2))')i,j,k,tHot(i,j,k),pgl(7,k,i,j)
!     *      ,pgl(3,k,i,j),cHot(i,j,k),proizv*razn,prion*(-pgi(6,k,i,j))
!	  end if
    3  continue
c*
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mass(23).ne.0) then
c    . . . Учет нелокального нагрева
         q_ot=pgl(19,n0,i,j)+pgl(19,n0,its-i+1,j)
         if(pgl(19,n0,i,j).ge.pgl(19,n0,its-i+1,j)) then
            q_max=pgl(19,n0,i,j)
         else
            q_max=pgl(19,n0,its-i+1,j)
         end if
         q_ot=q_ot/q_max
         pot_m=0.2
	 ! pot_m=0.02                    !   Поток эрг/см2*c-1
c  Задание ночного потока тепла в эрг/см2*с-1 на верхней границе
         if (klik.eq.0) then
           print*,' Поток тепла THot=',pot_m,' эрг/см2*c-1'
           klik=1
         end if
         pot=pot_m*.5*q_ot
         pot=pot/alyam(3)
         rp=rads(n0)-rads(n0-1)
         Thot(i,j,n0)=(pb(n0)+pot*rp)/(1.-pa(n0))
       else
         Thot(i,j,n0)=pb(n0)/(1.-pa(n0))
       end if
       do 4 l=2,n0m
        k=n0-l+1
        Thot(i,j,k)=pa(k+1)*Thot(i,j,k+1)+pb(k+1)
       !!!   if(Thot(i,j,k).lt.0) print*,i,j,k,Thot(i,j,k)

    4  continue
       ikey=0
       do k=1,nh
         if(Thot(i,j,k).lt.pgl(7,k,i,j))then
            Thot(i,j,k)=pgl(7,k,i,j)
	 end if
       end do
    2 continue
    1 continue
 ! 100 format('+ Thot=',f6.0,' k=',i3,' i=',i2,' j=',i2)
  100 format(' Thot=',f6.0,' k=',i3,' i=',i2,' j=',i2)
      
!      close (444)
      deallocate (pa,pb) 
      return

      end
C
!!!!  
   	subroutine thot_ic(an6,pgl,rads,kpars,n,n1,n2,dt,n0)
c . . .  . . . cyclic prog along meridian
c     . . .    nm=n1+n1-2
  
      dimension an6(n1,n2,n),pgl(kpars,n,n1,n2),rads(n)
      allocatable tn(:),v(:),a(:)
     *           ,b(:),c(:),f(:)
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
              v(i)=pgl(11,k,i,j)
              tn(i)=an6(i,j,k)
            end do
            do i=n1+1,nm
              v(i)=-pgl(11,k,i-n1+1,j+n_d)
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
      deallocate (tn,v,a,b,c,f)

      return
      end
      
	subroutine thot_jc(an6,pgl,rads,kpars,n,n1,n2,dt,n0)
c     . . . циклическая прогонка
    
       dimension an6(n1,n2,n),pgl(kpars,n,n1,n2),rads(n)
       allocatable tn(:),a(:),b(:),c(:),f(:)
       allocate (tn(n2),a(n2),b(n2),c(n2),f(n2))
       data pi/3.1415926/,re/6.371e8/,par/1./
       
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
            vm=(pgl(12,k,i,j)-abs(pgl(12,k,i,j)))*.5
            vp=(pgl(12,k,i,j)+abs(pgl(12,k,i,j)))*.5
            a(j)=vp*ot*par/del
            c(j)=-vm*ot*par/del
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
       deallocate (tn,a,b,c,f)
       return
       end
!!!!  
   	subroutine thot_icl(an6,cHot,pgl,rads,kpars,n,n1,n2,dt,n0)
c . . .  . . . cyclic prog along meridian
c . . .  . . . conductivity
c     . . .    nm=n1+n1-2
    
      
      dimension an6(n1,n2,n),cHot(n1,n2,n),pgl(kpars,n,n1,n2),rads(n)
     
      allocatable tn(:),v(:),alyam(:),rc(:)
     *,          a(:),b(:),c(:),f(:)
      data pi/3.1415926/,re/6.371e8/,par/1./ ,bk/1.38e-16/
      nm=n1+n1-2
      allocate (tn(nm),v(nm),alyam(nm),rc(nm)
     *,          a(nm),b(nm),c(nm),f(nm))
      ns=n1-1
      n_d=n2/2
      dtet=pi/ns
      ot=dt/dtet
      do k=2,n0
        rk=rads(k)+re
       do j=1,n_d
c . . . переход к меридиональному кругу
            do i=1,n1
              v(i)=pgl(11,k,i,j)
              tn(i)=an6(i,j,k)
	        rc(i)=1.5*cHot(i,j,k)*bk
	
c . . . conductivity
              alyam(i)=75.9*tn(i)**0.69
            end do
            do i=n1+1,nm
              v(i)=-pgl(11,k,i-n1+1,j+n_d)
              tn(i)=an6(i-n1+1,j+n_d,k)
	        rc(i)=1.5*cHot(i-n1+1,j+n_d,k)*bk
	
			alyam(i)=75.9*tn(i)**0.69
            end do
            do i=1,nm
	       
              im=i-1
              ip=i+1
              if(i.eq.nm) ip=1
              if(i.eq.1) im=nm
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
              a(i)=(vp*par+am_2/(del*rc(i)))*ot/rk
              c(i)=(-vm*par+ap_2/(del*rc(i)))*ot/rk
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

      deallocate (tn,v,alyam,rc,a,b,c,f)
      return
      end
      
	subroutine thot_jcl(an6,cHot,pgl,rads,kpars,n,n1,n2,dt,n0)
c     . . . cyclic prog along
c . . .   . conductivity
       
       dimension an6(n1,n2,n),cHot(n1,n2,n),pgl(kpars,n,n1,n2),rads(n)
       allocatable tn(:),alyam(:),a(:),b(:),c(:),f(:)
       data pi/3.1415926/,re/6.371e8/,par/1./,bk/1.38e-16/
       
       nm=n2
       allocate (tn(nm),alyam(nm),a(nm),b(nm),c(nm),f(nm))
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
c . . . conductivity
           if(tn(j).le.0) print*,tn(j),i,j,k
		 alyam(j)=75.9*tn(j)**0.69
          end do
          do j=1,n2
            rc=1.5*cHot(i,j,k)*bk
            jm=j-1
            jp=j+1
            if(j.eq.1) jm=n2
            if(j.eq.n2) jp=1
            vm=(pgl(12,k,i,j)-abs(pgl(12,k,i,j)))*.5
            vp=(pgl(12,k,i,j)+abs(pgl(12,k,i,j)))*.5
            ap_2=(alyam(jp)+alyam(j))*.5
            am_2=(alyam(j)+alyam(jm))*.5
            a(j)=(vp*par+am_2/(del*dfi*rc))/del*ot
            c(j)=(-vm*par+ap_2/(del*dfi*rc))/del*ot
            b(j)=1.+a(j)+c(j)
            f(j)=(tn(j)-vm*(1.-par)*ot/del*(tn(jp)-tn(j))-
     *                 vp*(1.-par)*ot/del*(tn(j)-tn(jm)))
          end do
          call cyclp(a,b,c,f,tn,n2)
          do j=1,n2
            an6(i,j,k)=tn(j)
          end do
        end do
       end do

       deallocate (tn,alyam,a,b,c,f)
       return
       end
