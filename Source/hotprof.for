!! EQUATION	
      subroutine hotprof(cHot,tHot,cNd,cNo,qHot,cNOi,cNe,pgl,pgi
     *           ,ctd,rads,rp,g,ins,kpars,nh,its,ids,i,j,hi,key,dt)
     
  	
      dimension pgl(kpars,nh,its,ids),pgi(ins,nh,its,ids),
     *          cHot(its,ids,nh),tHot(its,ids,nh),cNd(its,ids,nh)
     *         ,qHot(its,ids,nh)
     *         ,cNOi(nh),rads(nh),rp(nh),g(nh),ctd(nh),cNo(nh),cNe(nh)
      allocatable a(:),b(:),d(:),e(:),r(:),cmd(:)
     *          ,sdim(:),cH(:),cNv(:) ! N2v
	  allocate (a(nh),b(nh),d(nh),e(nh),r(nh),cmd(nh)
     *         ,sdim(nh),cH(nh),cNv(nh)) ! 
      INCLUDE 'alpha.inc'	 
      data amo2,amn2,amo,amno/ 53.12e-24,46.51e-24,26.56e-24,49.82e-24/
     *    ,bk/1.38e-16/
      data erg/1.6e-12/
      data ea1,      ea2,     eNv
     *   /0.18,     1.28,    0.19/,
     *     e7,       e9 ,     e12
     *   /2.45,     1.11,    3.58/ ,
     *     e17,       e19
     *   /0.68,      1.03/
      const=bk/amno
      const1=bk/amo2
      cH=cHot(i,j,1:nh)
    	do  k=2,nh-1
! . .  p - lost and...
         u=2.6e-09*pgl(3,k,i,j)*sqrt(0.79/8.)
!       	   
	   Gam=pgi(6,k,i,j)+pgl(7,k,i,j)
	  p=u+4.8e-13*sqrt(Gam)*pgi(1,k,i,j)*(10.6-0.67*alog10(Gam))**2
! . . . q - sourse O2 - Nd
         w=alyam7*pgl(1,k,i,j)
	  qHot(i,j,k)=w*e7 !heat energy in ev
! . . .  NO - Nd
	  tmp=alyam12*cNo(k)
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
	    !!! N2v
	    Tv=pgl(9,k,i,j)	!tv=te
	    pok=3353.0/Tv
	    ! if (pok.gt.60.) pok=60.
           est=exp(-pok)
         
          cNv(k)= pgl(2,k,i,j)*(1.0-est)*est
	    pok1=69.9/pgl(7,k,i,j)**0.333
	    tmp=cNv(k)*pgl(3,k,i,j)*1.07e-10*exp(-pok1)
          q=w+tmp
          qHot(i,j,k)=(qHot(i,j,k)+tmp*eNv)*erg ! energi in erg
! . . .        
       ! EQUATION
         pro= (rp(k)+rp(k-1))*.5
         a(k)=(1./dt+p)*pro
         b(k)=(q+cHot(i,j,k)/dt)*pro
   !     if(k.ne.nh) then
           dTh=tHot(i,j,k+1)-tHot(i,j,k-1)
           crab=0.5*dTh/(tHot(i,j,1)*pro)
   !      else
   !        dTh=tHot(i,j,k)-tHot(i,j,k-1)
   !        crab=dTh/(tHot(i,j,k)*rp(k))
    !     end if
c     . . . scale height
         h =const*tHot(i,j,k)/g(k)	 ! test 
         sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
         ams=(amo*pgl(3,k,i,j)+amo2*pgl(1,k,i,j)+
     *        amn2*pgl(2,k,i,j))/sum
         hsr=bk*pgl(7,k,i,j)/(ams*g(k))
c    . . . Coef. Mol. Dif.
         cmd(k)=3.e17/sum*sqrt(pgl(7,k,i,j)) ! test for Km
         vrp=abs(pgl(10,k+1,i,j))
         vrm=abs(pgl(10,k-1,i,j))
         r(k)=cmd(k)*(1./h+crab)+ctd(k)*(1./hsr+crab)
         r(k)=r(k)-pgl(10,k,i,j)
  !       r(k)=r(k)-(pgl(10,k+1,i,j)+vrp)*.5-(pgl(10,k-1,i,j)-vrm)*.5
        end do
        a(1)=(rp(1)/dt)
        b(1)=(cHot(i,j,1)/dt)*rp(1)
        sum=pgl(1,1,i,j)+pgl(2,1,i,j)+pgl(3,1,i,j)
        cmd(1)=3.e17/sum*sqrt(tHot(i,j,1))
        dTh=tHot(i,j,2)-tHot(i,j,1)
        crab=dTh/tHot(i,j,1)/rp(1)
        h =const*tHot(i,j,1)/g(1)
        ams=(amo*pgl(3,1,i,j)+amo2*pgl(1,1,i,j)+
     *        amn2*pgl(2,1,i,j))/sum
        hsr=bk*pgl(7,1,i,j)/(ams*g(1))
        r(1)=cmd(1)*(1./h+crab)+ctd(1)*(1./hsr+crab)
        do k=1,nh-1
         cmd05=(cmd(k+1)+cmd(k)+ctd(k+1)+ctd(k))*.5
         r05=(r(k)+r(k+1))*.5
  !       if(k.eq.nh-1) then
  !         h =const*pgl(7,nh-1,i,j)/g(nh-1)  !test
  !         r05=cmd(nh-1)/h
  !       end if
         d(k+1)=cmd05/rp(k)+r05*.5
         e(k+1)=(cmd05/rp(k)-r05*.5)/d(k+1)
         if(e(k+1).lt.0.) then
          d(k+1)=cmd05/rp(k)+r(k+1)
          e(k+1)=(cmd05/rp(k))/d(k+1)
         end if
      end do
      call potno(cH,sdim,a,b,d,e,nh) !test
      cHot(i,j,1:nh)=ch
	  deallocate (a,b,d,e,r,cmd
     *          ,sdim,cH,cNv) 
	  return
  	  end
