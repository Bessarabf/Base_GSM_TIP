      subroutine hotprog(cHot,tHot,cNd,qHot,cNOi,cNe,pgl,pgi
     *                 ,ctd,rads,rp,g,ins,kpars,nh,its,ids,i,j,dt)
     
      dimension pgl(kpars,nh,its,ids),pgi(ins,nh,its,ids),
     *          cHot(its,ids,nh),tHot(its,ids,nh),cNd(its,ids,nh)
     *         ,qHot(its,ids,nh)
     *         ,cNOi(nh),cNe(nh),ctd(nh),RADS(NH),rp(nh),g(nh)
      allocatable a(:),b(:),c(:),f(:),cmd(:)
     *         ,h(:),hsr(:),alf(:),bet(:)
     *         ,cNv(:),cNv2(:) ! N2v
     *         ,pa(:),pb(:) 
      allocate (a(nh),b(nh),c(nh),f(nh),cmd(nh)
     *         ,h(nh),hsr(nh),alf(nh),bet(nh)
     *         ,cNv(nh),cNv2(nh)                 ! N2v
     *         ,pa(nh+1),pb(nh+1)) 
      INCLUDE 'alpha.inc'	 
      data amo2,amn2,amo/ 53.12e-24,46.51e-24,26.56e-24/
     *    ,bk/1.38e-16/,re/6.371e8/,pi/3.14159/
      data erg/1.6e-12/
      data ea1,      ea2,     eNv
     *   /0.18,     1.28,    0.19/,
     *     e7,       e9 ,     e12
     *   /2.45,     1.11,    3.58/ ,
     *     e17,       e19
     *   /0.68,      1.03/
      
      const=bk/amo
      dtet=pi/(its-1)
      dfi=2.*pi/ids
      teta=dtet*(i-1)
      sin_t=sin(teta)
      cot_t=cos(teta)/sin_t
      jp=j+1
      jm=j-1
      if(j.eq.ids) jp=1
      if(j.eq.1)  jm=ids
      
	pa(2)=0
	pb(2)=cHot(i,j,1)
      do k=1,nh  
          
c     . . . scale height
         h(k) =const*Thot(i,j,k)/g(k)
         sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
         ams=(amo*pgl(3,k,i,j)+amo2*pgl(1,k,i,j)+
     *        amn2*pgl(2,k,i,j))/sum
         hsr(k)=bk*pgl(7,k,i,j)/(ams*g(k))
c    . . . Coef. Mol. Dif.
         cmd(k)=1.85e18/sum*sqrt(thot(i,j,k)) ! 10.03.2010
!        cmd(k)=4.55e17/sum*sqrt(thot(i,j,k)) ! test for Km

!	   cmd(k)=4.55e17/sum*sqrt(PGL(7,K,I,J)) ! test for Km

         alf(k)=cmd(k)/(cmd(k)+ctd(k))
         bet(k)=ctd(k)/(cmd(k)+ctd(k))
      end do	
      do  k=2,nh-1
            
         rk=rads(k)+re
         pro= (rp(k)+rp(k-1))*.5
! . .  p - lost and...
!        u=2.6e-09*pgl(3,k,i,j)*sqrt(0.79/8.)
         u=1.143e-12*(pgl(3,k,i,j)*sqrt(thot(i,j,k)+pgl(7,k,i,j))+!hard sphere 
     *   pgl(2,k,i,j)*2.031*sqrt(thot(i,j,k)+pgl(7,k,i,j)*0.571))!approx:Brunel  
!     for Oh - O+  Namgaladze 
  	   Gam=pgi(6,k,i,j)+pgl(7,k,i,j)
!!	    p=u+4.8e-13*sqrt(Gam)*pgi(1,k,i,j)*(10.6-0.67*alog10(Gam))**2
!                  Stubbe
         p=u+1.86e-9*pgi(1,k,i,j)*(Gam/2000.)**0.37
! . . .  q - sourse O2 - Nd
        w=alyam7*pgl(1,k,i,j)
        qHot(i,j,k)=w*e7 !heat energy in ev
! . . .  NO - Nd
      
	  tmp=alyam12*pgl(4,k,i,j)
         w=w+tmp
         qHot(i,j,k)=qHot(i,j,k)+tmp*e12 !heat energy in eV
! . . .  O - Nd
         tmp=alyam9*pgl(3,k,i,j)
         w=w+tmp
         qHot(i,j,k)=qHot(i,j,k)+tmp*e9
! . . . 
         tmp=alyam17*pgi(1,k,i,j)
         w=w+tmp
         qHot(i,j,k)=qHot(i,j,k)+tmp*e17
! . . . 
	     w=w*cNd(i,j,k)
         qHot(i,j,k)=qHot(i,j,k)*cNd(i,j,k)
	!!! O2 - O+
         tmp=2.1e-11*pgl(1,k,i,j)*pgi(1,k,i,j)
	     w=w+tmp
         qHot(i,j,k)=qHot(i,j,k)+tmp*e19
        !!! NO+ - e
         otn=exp(300/pgl(9,k,i,j))
	   tmp=alfa1*otn*cNOi(k)*cNe(k)
         w=w+tmp 
         qHot(i,j,k)=qHot(i,j,k)+tmp*(r1*ea1+(1-r1)*ea2)
      !!! N2v
  	    Tv=1.3*pgl(7,k,i,j)	!tv=1.2tn
      ! 	Tv=1.15*pgl(9,k,i,j)	!tv=1.15te

	    pok=3353.0/Tv
	    ! if (pok.gt.60.) pok=60.
           est=exp(-pok)
      ! first vibrational laier    
          cNv(k)= pgl(2,k,i,j)*(1.0-est)*est
	    pok1=69.9/pgl(7,k,i,j)**0.333
	    tmp=cNv(k)*pgl(3,k,i,j)*1.07e-10*exp(-pok1)
          w=w+tmp
	    qHot(i,j,k)=(qHot(i,j,k)+tmp*eNv) 
	! second vibrational laier
	    pok2=2.*pok
	    cNv2(k)=pgl(2,k,i,j)*(1.0-est)*exp(-pok2)
	    tmp=cNv2(k)*pgl(3,k,i,j)*1.07e-10*exp(-pok1)
	    q=w+tmp
            qHot(i,j,k)=(qHot(i,j,k)+2*tmp*eNv)*erg
! . . .        
! EQUATION
      
!          dtnp=alog(pgl(7,k+1,i,j)/pgl(7,k,i,j))
!          dtnm=alog(pgl(7,k,i,j)/pgl(7,k-1,i,j))
           dtNp=alog(Thot(i,j,k+1)/Thot(i,j,k))
           dtNm=alog(Thot(i,j,k)/Thot(i,j,k-1))
c      . . . Є-в ¤Ёддг§ЁЁ ў ¤а®Ў­®© в®зЄҐ
          cmdp=(cmd(k+1)+cmd(k)+ctd(k+1)+ctd(k))*.5
          cmdm=(cmd(k)+cmd(k-1)+ctd(k)+ctd(k-1))*.5

          a(k)=(cmdp/pro)/rp(k)
          c(k)=(cmdm/pro)/rp(k-1)
          b(k)=a(k)+c(k)+1./dt+p
	
          clp= 0.5*(alf(k+1)/h(k+1)+dtnp/rp(k)+bet(k+1)/hsr(k+1))  !!! thot=tn
          clm=-0.5*(alf(k-1)/h(k-1)+dtnm/rp(k-1)+bet(k-1)/hsr(k-1))!!!!thot=tn
          aprim=(cmd(k+1)+ctd(k+1))/pro*clp
          cprim=(cmd(k-1)+ctd(k-1))/pro*clm
          cl=0.5*(-dtnp/rp(k)+dtnm/rp(k-1))
          a(k)=a(k)+aprim -pgl(10,k+1,i,j)/pro*.5 !   
          c(k)=c(k)+cprim +pgl(10,k-1,i,j)/pro*.5
          bprim=cl*(cmd(k)+ctd(k))/pro
		b(k)=b(k)+bprim !!!!
      
		f(k)=q+cHot(i,j,k)/dt
c . . .  „ЁўҐаЈҐ­жЁп V:
          del=2.*dfi*rk*sin_t
          div=(pgl(12,k,i,jp)-pgl(12,k,i,jm))/del
c . . .  гзҐв Є®в ­ЈҐ­б :
          del=2.*dtet*rk
          div=div+pgl(11,k,i,j)/rk*cot_t+
     *       (pgl(11,k,i+1,j)-pgl(11,k,i-1,j))/del
      	if(div.gt.0.) then
             b(k)=b(k)+div
          else
             f(k)=f(k)-div*cHot(i,j,k)
          end if
	! . . .  Прямая прогонка
          del1= (b(k)-pa(k)*c(k))
		pa(k+1)=a(k)/del1
!		if(pa(nh).gt.1.) then
!	print*,pa(nh),pb(nh),aprim,bprim,cprim,b(k),b(k-1),div
!	pause
!	end if

          pb(k+1)=(c(k)*pb(k)+f(k))/del1
	

        end do
c . . . элемент массива используется для удобства
c       
!        f(nh)=exp(-rp(nh-1)/h(nh))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Gershengorn
!
	  f(nh)=(2*h(nh)-rp(nh-1))/(rp(nh-1)+2*h(nh))
	  pot_m=0.006 ! 0.02 from cmd(nh)
	  !/cmd(nh)
	  !pot=0.
	       q_ot=pgl(18,nh,i,j)+pgl(18,nh,its-i+1,j)
         if(pgl(18,nh,i,j).ge.pgl(18,nh,its-i+1,j)) then
            q_max=pgl(18,nh,i,j)
         else
            q_max=pgl(18,nh,its-i+1,j)
         end if
         q_ot=q_ot/q_max
	   pot=pot_m*.5*q_ot
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !! pot=0.							 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ger=(2*h(nh)*rp(nh-1))/(rp(nh-1)+2*h(nh))

        cHot(i,j,nh)=(f(nh)*pb(nh)+pot*ger)/(1.-pa(nh)*f(nh))
        
		   
	! обратная прогонка
	do k=nh-1,2,-1
          cHot(i,j,k)=pa(k+1)*cHot(i,j,k+1)+pb(k+1)
		
	    if(cHot(i,j,k).le.1.) cHot(i,j,k)=1.
	end do
      
	deallocate (a,b,c,f,cmd
     *         ,h,hsr,alf,bet,cNv,cNv2,pa,pb)


      return

  	end
