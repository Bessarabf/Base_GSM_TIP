c. . . programm list
c - nonew    1
c - nnew     2
c - nnewp    3
c - ndnew    4
c - conn2i   5
c - cono2i   6
c - connoi   7
c - chep     8
c - potno    9
c - hori    10
c - horj    11
 
c
      subroutine nonew(pgl,pgi,gkoor,ctd,rads,rp,g,
     *                 kpars,ins,nh,its,ids,del,uts,dt,mass)
!!!!     
        USE mo_bas_gsm
!!!!
	dimension pgl(kpars,nh,its,ids),pgi(ins,nh,its,ids),
     *          mass(30),

     *          rads(nh),rp(nh),g(nh),ctd(nh),gkoor(2,its,ids)
     
     	
	allocatable cNd(:,:,:),cO2i(:),cNoi(:)
     *           ,cN2i(:),cNe(:),cNo(:)
	allocate (cNd(its,ids,nh),cO2i(nh),cNoi(nh)
     *           ,cN2i(nh),cNe(nh),cNo(nh))
	           
      INCLUDE 'alpha.inc'

      data key/0/
      cr=180.d0/pi
      bkg=bk/amo2/g(1)
      do i=1,its
       do j=1,ids 
 	   do k=1,nh
         cNd(i,j,k)=pgl(17,k,i,j)   

	   end do
	 end do
	end do
	do i=2,its-1
       do j=1,ids
c      . . . calculation of cos(hi)
        rlat=gkoor(1,i,j)/cr
        rlat=pi/2.-rlat
        dolgg=gkoor(2,i,j)/cr
        cx=sin(del)*sin(rlat)+cos(del)*cos(rlat)
     *    *cos(om*(uts-43200.)+dolgg)
        hi=acos(cx)
        do k=1,nh
         cNe(k)=pgl(6,k,i,j)+pgi(1,k,i,j)
         cNo(k)=pgl(4,k,i,j)
        end do
        !!! smoothing sudden change Ne
	  cNe(16)=0.5*(cNe(15)+cNe(17))
!!!!!!!!!!! photochem approx
	  
        if ((mass(20).ne.2)) then
           !!! photochem approx
           call cono2i (cO2i,cNo,cNe,pgl,kpars,nh,its,ids,i,j)
           call connoi (cNoi,cO2i,cN2i,cNe,pgl,rp,g,
     *                  kpars,nh,its,ids,i,j)
           call conn2i (cN2i,cNe,pgl,kpars,nh,its,ids,i,j)
           do k=1,nh    
            pgl(18,k,i,j)=cO2i(k)
            pgl(19,k,i,j)=cNOi(k)
           
           end do   
        else
           do k=1,nh    
            cO2i(k)=pgl(18,k,i,j)
            cNOi(k)=pgl(19,k,i,j)
            !cN2i(k)=pgl(6,k,i,j)-cNOi(k)-cO2i(k)
	      call conn2i (cN2i,cNe,pgl,kpars,nh,its,ids,i,j)
           end do     
        end if

46      format(1p4e8.1)
! . . . cNd - calculation
!  photochen N2d       
!        call ndnew (cNd,cNo,cNoi,cN2i,cNe
!     *             ,pgl,kpars,nh,its,ids,i,j,key,dt)
! 
!  progonka N(2D)
          call ndprog(cNd,cNo,cNoi,cN2i,cNe
     *            ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
!  progonka N(4S)
       
    	  call nprog(cNd,cNoi,cO2i,cNe
     *            ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
! progonka NO. Altitude part
	   call noprog(cNd,cNo,cO2i
     *            ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
          end do
	end do	  
      !  S-N poles smoothing for N2D
	call bongl(cNd,nh,its,ids)
      !  horisontal part for NO
        call bospgl(pgl,kpars,nh,its,ids,4)
	call pgl_ic(pgl,rads,kpars,nh,its,ids,dt,nh,4)
        call pgl_jc(pgl,rads,kpars,nh,its,ids,dt,nh,4) 
        call bospgl(pgl,kpars,nh,its,ids,4)
      !  horisontal part for N(4s)
	call bospgl(pgl,kpars,nh,its,ids,5)
	call pgl_ic(pgl,rads,kpars,nh,its,ids,dt,nh,5)
        call pgl_jc(pgl,rads,kpars,nh,its,ids,dt,nh,5)
	call bospgl(pgl,kpars,nh,its,ids,5)
      !  horisontal part for N(2D)
      	DO K=1,NH
	   DO I=1,ITS
	      DO J=1,IDS
	       PGL(17,K,I,J)=CND(I,J,K)
            
              END DO
	   END DO
	END DO
       	call pgl_ic(pgl,rads,kpars,nh,its,ids,dt,nh,17)
        call pgl_jc(pgl,rads,kpars,nh,its,ids,dt,nh,17)
	call bospgl(pgl,kpars,nh,its,ids,17)

c     . . .  horisontal circulation NO,N
!      call hori(pgl,rads,kpars,nh,its,ids,4,dt)
!      call horj(pgl,rads,kpars,nh,its,ids,4,dt)
!      call hori(pgl,rads,kpars,nh,its,ids,5,dt)
!      call horj(pgl,rads,kpars,nh,its,ids,5,dt)
!
! if photochem approx 
      if ((mass(20).ne.2)) then
          call bospgl(pgl,kpars,nh,its,ids,18)
          call bospgl(pgl,kpars,nh,its,ids,19)
      end if
      deallocate (cNd,cO2i,cNoi,cN2i,cNe,cNo)
      return
      end

c
c . . . for GSM [N4S)] диффузионная постановка
      subroutine nnewp(cNd,cNo,cNoi,cO2i,cNe
     *                ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
      dimension pgl(kpars,nh,its,ids),cNd(its,ids,nh),
     *          cNoi(nh), cO2i(nh),cNe(nh),cNo(nh)
     *         ,ctd(nh),rads(nh),rp(nh),g(nh)
     *         
      allocatable a(:),b(:),d(:),e(:),r(:),cmd(:)
     *         ,sdim(:),cN(:)
      allocate(a(nh),b(nh),d(nh),e(nh),r(nh),cmd(nh)
     *         ,sdim(nh),cN(nh))
      INCLUDE 'alpha.inc'
      data  am o2,amn2,amo,amn/ 53.12e-24,46.51e-24,26.56e-24
     *                         ,23.26e-24/
     *     ,bk/1.38e-16/, re/ 6.371e+8/
c
c     data  alfa1/4.2e-7/,r1,r3/0.7,0.5/
c    *    , alyam2,  alyam8,   alyam9
c    *    /1.2e-10, 4.4e-12,  4.5e-13/
c    *    ,alyam11, alyam10,  alyam13, alyam14, alyam16
c    *    /3.6e-10, 3.4e-11,  1.06e-5, 4.5e-6,  2.3e-14/
c    *    ,alyam6 / 0.1/
c
      const=bk/amn
      const1=bk/amo2
      cN(1)=pgl(5,1,i,j)
      do 1 k=2,nh-1
        ot=300./pgl(9,k,i,j)
        ots=sqrt(ot)
c
c p -lost, and q - source
        u=alfa1*ot**0.85*(1.-r1)*cNoi(k)*cNe(k)
        u=u+(alyam9*pgl(3,k,i,j)+alyam11/ots*cNe(k)+
     *    alyam16*pgl(2,k,i,j)+alyam13)*cNd(i,j,k)
        q=u+r3*alyam6*pgl(14,k,i,j)
        h o2=const1*pgl(7,k,i,j)/g(k)
        x=(rads(k)+re)/h o2
        ch=chep(x,hi)
        pok=1.e-8*(pgl(1,k,i,j)*h o2*ch*1.E-4)**0.38
        if (pok.gt.60.) pok=60.
        q=q+alyam14*exp(-pok)*cNo(k)
c . . .  lost
        w=alyam2*cO2i(k)
        w=w+alyam8*exp(-3220./pgl(7,k,i,j))*pgl(1,k,i,j)
        p=w+alyam10*cNo(k)
        pro= (rp(k)+rp(k-1))*.5
        a(k)=(1./dt+p)*pro
        b(k)=(q+pgl(5,k,i,j)/dt)*pro
        dtn=pgl(7,k,i,j)-pgl(7,k-1,i,j)
        crab=dtn/pgl(7,k,i,j)/rp(k)
c     . . . scale height
        h =const*pgl(7,k,i,j)/g(k)
        sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
        ams=(amo*pgl(3,k,i,j)+amo2*pgl(1,k,i,j)+
     *       amn2*pgl(2,k,i,j))/sum
        hsr=bk*pgl(7,k,i,j)/(ams*g(k))
c    . . . Coef. Mol. Dif.
        cmd(k)=4.55e17/sum*sqrt(pgl(7,k,i,j))
        vrp=abs(pgl(10,k,i,j))
        vrm=abs(pgl(10,k-1,i,j))
        r(k)=cmd(k)*(1./h+crab)+ctd(k)*(1./hsr+crab)
       
	  r(k)=r(k)-(pgl(10,k+1,i,j)+vrp)*.5-(pgl(10,k-1,i,j)-vrm)*.5
cc      r(k)=r(k)-pgl(10,k,i,j)
   1  continue
      a(1)=(rp(1)/dt)
      b(1)=(cN(1)/dt)*rp(1)
      sum=pgl(1,1,i,j)+pgl(2,1,i,j)+pgl(3,1,i,j)
      cmd(1)=4.55e17/sum*sqrt(pgl(7,1,i,j))
      dtn=pgl(7,2,i,j)-pgl(7,1,i,j)
      crab=dtn/pgl(7,1,i,j)/rp(1)
      h =const*pgl(7,1,i,j)/g(1)
      ams=(amo*pgl(3,1,i,j)+amo2*pgl(1,1,i,j)+
     *     amn2*pgl(2,1,i,j))/sum
      hsr=bk*pgl(7,1,i,j)/(ams*g(1))
      r(1)=cmd(1)*(1./h+crab)+ctd(1)*(1./hsr+crab)
      do k=1,nh-1
        cmd05=(cmd(k+1)+cmd(k)+ctd(k+1)+ctd(k))*.5
        r05=(r(k)+r(k+1))*.5
c       if(k.eq.n-1) then
c        h =const*pgl(7,n-1,i,j)/g(n-1)
c        r05=cmd(n-1)/h
c      end if
        d(k+1)=cmd05/rp(k)+r05*.5
        e(k+1)=(cmd05/rp(k)-r05*.5)/d(k+1)
        if(e(k+1).lt.0.) then
           d(k+1)=cmd05/rp(k)+r(k+1)
           e(k+1)=(cmd05/rp(k))/d(k+1)
        end if
       end do
       call potno(cN,sdim,a,b,d,e,nh)
       do k=2,nh
         pgl(5,k,i,j)=cN(k)
       end do
	deallocate(a,b,d,e,r,cmd,sdim,cN)

      return
      end
c

c . . . for GSM
      subroutine ndnew(cNd,cNo,cNoi,cN2i,cNe
     *                ,pgl,kpars,nh,its,ids,i,j,key,dt)
      dimension cNd(its,ids,nh),cNoi(nh),cN2i(nh),cNe(nh),
     *          cNo(nh),pgl(kpars,nh,its,ids)
      INCLUDE 'alpha.inc'
c
c     data  alfa1/4.2e-7/,r1,r3/0.7,0.5/
c    *    , alyam4 /1.4e-10/
c    *    , alyam6,  alyam7,   alyam9
c    *    /    0.1,  5.e-12,  4.5e-13/
c    *    ,alyam11, alyam12,  alyam13, alyam16
c    *    /3.6e-10,  7.e-11,  1.06e-5, 2.3e-14/

      do 1 k=1,nh
       ot=300./pgl(9,k,i,j)
       ots=sqrt(ot)
       tr=(pgl(8,k,i,j)+pgl(7,k,i,j))*.5
c p -lost, and q- source
c
       a=alfa1*ot**0.85*r1*cNoi(k)*cNe(k)
       a=a+alyam4*(300./tr)**0.44*cN2i(k)*pgl(3,k,i,j)
       q=a+r3*alyam6*pgl(14,k,i,j)
       b=alyam7*pgl(1,k,i,j)
       b=b+alyam9*pgl(3,k,i,j)
       b=b+alyam11/ots*cNe(k)
       b=b+alyam12*cNo(k)+alyam16*pgl(2,k,i,j)
       p=b+alyam13
      ! if(key.eq.0) then
      !  cNd(i,j,k)=q/p
      ! else
        cNd(i,j,k)=(cNd(i,j,k)+q*dt)/(1.+p*dt)
      ! end if
   1  continue
      return
      end
c
c . . . for GSM
      subroutine nnew(cNd,cNo,cNoi,cO2i,cNe
     *                ,pgl,rads,g,kpars,nh,its,ids,i,j,hi,dt)
      dimension pgl(kpars,nh,its,ids),cNd(its,ids,nh),
     *          cNoi(nh), cO2i(nh),cNe(nh),cNo(nh),rads(nh),g(nh)
      data  am o2,amn2,amo,amn/ 53.12e-24,46.51e-24,26.56e-24
     *                         ,23.26e-24/
     *     ,bk/1.38e-16/, re/ 6.371e+8/
      data  alfa1/4.2e-7/,r1,r3/0.7,0.5/
     *    , alyam2,  alyam8,   alyam9
     *    /1.2e-10, 4.4e-12,  4.5e-13/
     *    ,alyam11, alyam10,  alyam13, alyam14, alyam16
     *    /3.6e-10, 3.4e-11,  1.06e-5, 4.5e-6,  2.3e-14/
     *    ,alyam6 / 0.1/
      const=bk/amn
      const1=bk/amo2
      do 1 k=1,nh
       ot=300./pgl(9,k,i,j)
       ots=sqrt(ot)
       tr=(pgl(8,k,i,j)+pgl(7,k,i,j))*.5
c
c p -lost, and q- sourсe
       a=alfa1*ot**0.85*(1.-r1)*cNoi(k)*cNe(k)
       a=a+(alyam9*pgl(3,k,i,j)+alyam11/ots*cNe(k)+
     *   alyam16*pgl(2,k,i,j)+alyam13)*cNd(i,j,k)
       q=a+r3*alyam6*pgl(14,k,i,j)
       h o2=const1*pgl(7,k,i,j)/g(k)
       x=(rads(k)+re)/h o2
       ch=chep(x,hi)
       pok=1.e-8*(pgl(1,k,i,j)*h o2*ch*1.E-4)**0.38
       if (pok.gt.60.) pok=60.
       q=q+alyam14*exp(-pok)*cNo(k)
       w=alyam2*cO2i(k)
       b=alyam2*cO2i(k)
       b=b+alyam8*exp(-3220./pgl(7,k,i,j))*pgl(1,k,i,j)
       p=b+alyam10*cNo(k)
       if (k.eq.1) then
        pgl(5,1,i,j)=5.e5
       else
        pgl(5,k,i,j)=(pgl(5,k,i,j)+q*dt)/(1.+p*dt)
       end if
   1  continue
      return
      end
c
c     . . . horizontal transport for NO, N
      subroutine horj(pgl,r,kpars,nh,its,ids,np,dt)
      dimension pgl(kpars,nh,its,ids)
     *         ,r(nh)
      allocatable as(:)
      allocate (as(ids))
      data pi/3.14159265359d0/,re/6.371e+8/
      
      it= its-1
      dfi=2.*pi/ids
      dtet=pi/it
      do 101 k = 2 , nh
        rk=r(k)+re
        do 10 i = 2, it
          tet=dtet*(i-1)
          sint=sin(tet)
          do 1 j = 1, ids
            jp=j+1
            jm=j-1
            if(j.eq.1) jm=ids
            if(j.eq.ids) jp=1
            s=alog(pgl(np,k,i,j))
            sp1=alog(pgl(np,k,i,jp))
            sm1=alog(pgl(np,k,i,jm))
            vjm= (pgl(12,k,i,j)-abs(pgl(12,k,i,j)))*.5
            vjpl=(pgl(12,k,i,j)+abs(pgl(12,k,i,j)))*.5
            grad=vjpl*(s-sm1)+vjm*(sp1-s)
            div = (pgl(12,k,i,jp)-pgl(12,k,i,jm))*.5
            f=(div+grad)/(rk*sint*dfi)
            as(j)=s-f*dt
    1     continue
          do j = 1 , ids
           pgl(np,k,i,j)=exp(as(j))
          end do
  10    continue
 101  continue
      deallocate (as)
      return
      end
c     . . . horizontal transport for NO, N
      subroutine hori(pgl,r,kpars,nh,its,ids,np,dt)

	dimension pgl(kpars,nh,its,ids)
     *         ,r(nh)
      allocatable as(:)
	allocate (as(its))
      data pi/3.14159265359d0/,re/6.371e+8/
      it= its-1
      dtet=pi/it
      do 101 k = 2 , nh
        rk=r(k)+re
        do 10 j = 1 , ids
          do 1 i = 2 , it
            tet=dtet*(i-1)
            cot=1./tan(tet)
            if(pgl(np,k,i+1,j).le.0) then
	        print*,i+1,j,k,np, pgl(np,k,i+1,j)
	      end if
		  s=alog(pgl(np,k,i,j))
            sp1=alog(pgl(np,k,i+1,j))
            sm1=alog(pgl(np,k,i-1,j))
            vim= (pgl(11,k,i,j)-abs(pgl(11,k,i,j)))*.5
            vipl=(pgl(11,k,i,j)+abs(pgl(11,k,i,j)))*.5
            grad=vipl*(s-sm1)+vim*(sp1-s)
            grad=grad/dtet
            div = (pgl(11,k,i+1,j)-pgl(11,k,i-1,j))*.5
            div = div/dtet+pgl(11,k,i,j)*cot
            f=(div+grad)/rk
            as(i)=s-f*dt
    1     continue
          do i = 2 , it
           pgl(np,k,i,j)=exp(as(i))
          end do
  10    continue
 101  continue
      deallocate (as)
      return
      end
c
      function chep(x,hi)
c     . . . CHEPMEN's function
c     definition f-formula
      data pi/3.14159/
      ch(a,b,hi)=a*(1./cos(b*hi)-0.834)
      hig=hi/pi*180.
      a=3.88/x**1.143
      sq=sqrt(x)
      b=1.0123-1.454/sq
      if(hig.lt.89.) then
       d=ch(a,b,hi)
       chep=1./cos(hi-d)
       return
      else
       z=sin(hi)*x
       sqz=sqrt(z)
       pok=x-z
       if(pok.gt.30.)pok=30.
       ex=exp(pok)
       prom=2.507*sqz+.903/sqz
       hip=pi-hi
       d=ch(a,b,hip)
       chep=ex*prom-1./cos(hip-d)
      endif
      return
      end
c
c . . . for GSM N2+
      subroutine conn2i(cN2i,cNe,pgl,kpars,nh,its,ids,i,j)
      dimension cN2i(nh),cNe(nh),pgl(kpars,nh,its,ids)
      data   alfa3
     *    /1.8e- 7/
     *    , alyam4,  alyam5
     *    /1.4e-10, 5.0e-11/
        do k=1,nh
         tr=(pgl(8,k,i,j)+pgl(7,k,i,j))*.5
         r=300./tr
         amu=alyam4*r**0.44*pgl(3,k,i,j)
         amu=amu+alyam5*r*pgl(1,k,i,j)
         del=alfa3*(300./pgl(9,k,i,j))**0.4*cNe(k)+amu
         cN2i(k)=pgl(14,k,i,j)/del
         if(cN2i(k).lt.1.e-2) cN2i(k)=1.e-2
        end do
      return
      end
c . . . for GSM NO+
      subroutine connoi(cNoi,cO2i,cN2i,cNe,pgl,rp,g,
     *                  kpars,nh,its,ids,i,j)
      dimension cO2i(nh),cNoi(nh),cN2i(nh),pgl(kpars,nh,its,ids)
     *          ,cNe(nh),rp(nh),g(nh)  
      data amno/49.82e-24/
	alfa1=4.2e-7
***************************************************************
  !!     do 1 k=1,nh
  !!!!       coef=0.8*alfa1*(300./pgl(9,k,i,j))**0.85*cNe(k)
c        cmoli=pgl(13,k,i,j)+pgl(14,k,i,j)+pgl(16,k,i,j)
c        cmoli=sqrt(cmoli/coef)
 !!!!        cmoli=pgl(6,k,i,j)
*         cNoi(k)=cmoli-cO2i(k)-cN2i(k)
*  !!!!       cNOi(k)=pgl(15,k,i,j)/coef
*  !!! 	  if(cNoi(k).le.0.) then
*          print 100,cNoi(k),i,j,k
*100       format('+NO+ < 0 in connoi =',1pe9.2,3i4)
c         cNoi(k)=abs(cNoi(k))
*          cNoi(k)=1.
*        end if
*    1 continue
	!!! more thin condition	!!!!!!!!!
	do k =1,nh
	  sum= cO2i(k)+cN2i(k)
	  cNOi(k)= pgl(6,k,i,j)-sum
	
        if (cNOi(k).le.(0.1*sum))then
	!!!!!!!   call baromi(cNOi,pgl,rp,g,amno,kpars,nh,its,ids,i,j,k)
	!!!!!!!   exit
	      cNOi(k)=0.1*sum
	  end if
	end do
	 !!! smoothing sudden change NOi
 	cNOi(16)=0.5*(cNOi(15)+cNOi(17))
	return
      end


c . . . for GSM O2+
      subroutine cono2i(cO2i,cNo,cNe,pgl,kpars,nh,its,ids,i,j)
      dimension cO2i(nh),cNo(nh),cNe(nh),pgl(kpars,nh,its,ids)
	
      data alyam1, alyam2,  alyam3
     *    /4.4e-10,1.2e-10, 5.0e-16/
      do k=1,nh
         amu=alyam1*cNo(k)
         amu=amu+alyam2*pgl(5,k,i,j)+alyam3*pgl(2,k,i,j)
         r=300./pgl(9,k,i,j)
         if(pgl(9,k,i,j).le.1200.) then
          r=r**0.85
          alfa1=2.7e-7*r
         else
          r=r**0.55
          alfa1=1.6e-7*r
         end if
         del=amu+alfa1*cNe(k)
         cO2i(k)=pgl(13, k,i,j)/del
         if(cO2i(k).lt.1.e-2) co2i(k)=1.e-2
      end do
      return
      end
c     . . . потоковая прогонка
      subroutine potno(dim,sdim,a,b,d,e,n)
c potokowa progonka
c     dim -
c     sdim -
c     a, b, d - koefficients
c     bet,gam- progon koefficients
c
      dimension dim(n),a(n),b(n),d(n),e(n),sdim(n)
     
	  allocatable bet(:),gam(:)
	  allocate (bet(n+1),gam(n+1))
c      chas=d(2)+a(2)
       chas=d(2)+a(1)
       bet(1)=-1./chas
       gam(1)=dim(1)*d(2)/chas
      n1=n-1
c  . . . direct run
      do 1 k=1,n1
       if (abs(d(k+1)).lt.1.) then
        c=1./ d(k+1)
        del      = 1.+e(k+1)*(c-bet(k))*a(k+1)
        bet(k+1) = -e(k+1)*(c-bet(k))/del
        gam(k+1) = e(k+1)*(gam(k)-bet(k)*b(k))/del
        if(abs(del).lt.1) then
         print 100,del,k
100     format('alf.gt.1',e8.1,i4)
        end if
       else
        if(abs(e(k+1)).gt.1.) then
          del= d(k+1)+e(k+1)*(1.-bet(k)*d(k+1))*a(k+1)
          if(abs(del).lt.abs(d(k+1))) then
           print *,' alf great 1',del,d(k+1),k
          end if
          gam(k+1)= e(k+1)*d(k+1)*(gam(k)-bet(k)*b(k))/del
          bet(k+1)= -e(k+1)*(1.-bet(k)*d(k+1))/del
        else
          f=1./e(k+1)
          del= d(k+1)*f+(1.-bet(k)*d(k+1))*a(k+1)
          if(abs(del).lt.abs(d(k+1))) then
           print *,' alf great 1',del,d(k+1),f,k
          end if
          gam(k+1)= d(k+1)*(gam(k)-bet(k)*b(k))/del
          bet(k+1)= -(1.-bet(k)*d(k+1))/del
        end if
       end if

  1   continue
      dim(n)=gam(n)/(1.+bet(n)*a(n))
!     print *,dim(n),gam(n),bet(n),a(n),b(n),d(n),e(n)
c     dim(n)= bet(n)*(-b(n)+gam(n))
c     dim(n)= bet(n)*(+b(n)-gam(n))
      do 2 k=n,2,-1
        if (k.eq.n) then
         sdim(k)=(1.+bet(k)*a(k))*b(k)-a(k)*gam(k)
        else
         sdim(k)=(1.+bet(k)*a(k))*(sdim(k+1)+b(k))-a(k)*gam(k)
         dim(k)=gam(k)-bet(k)*(sdim(k+1)+b(k))
		end if
!	  write(*,*) k,dim(k)

    2 continue
	  deallocate(bet,gam)

      return
      end 

	! NO profile
      subroutine noprof(cNd,cNo,cO2i
     *             ,pgl,ctd,rads,rp,g,kpars,nh,its,ids,i,j,hi,dt)
     
      dimension pgl(kpars,nh,its,ids),
     *          cO2i(nh),cNo(nh),cNd(its,ids,nh),
     *          rads(nh),rp(nh),g(nh),ctd(nh)
      allocatable a(:),b(:),d(:),e(:),r(:),cmd(:)
     *          ,sdim(:)

      INCLUDE 'alpha.inc'	 
      data amo2,amn2,amo,amno/ 53.12e-24,46.51e-24,26.56e-24,49.82e-24/
     *    ,bk/1.38e-16/
      data re/6.371e+8/
	allocate (a(nh),b(nh),d(nh),e(nh),r(nh),cmd(nh)
     *          ,sdim(nh))

      const=bk/amno
      const1=bk/amo2
        do 1 k=2,nh-1
! . .  q - sourse and...
         u=(alyam7*cNd(i,j,k)+
     *      alyam8*exp(-3220./pgl(7,k,i,j))*pgl(5,k,i,j))*pgl(1,k,i,j)
         q=u+alyam3*cO2i(k)*pgl(2,k,i,j)
! . . . p - lost
         w=alyam10*pgl(5,k,i,j)+alyam12*cNd(i,j,k)
         w=w+alyam15*pgl(3,k,i,j)*pgl(2,k,i,j)*exp(940./pgl(7,k,i,j))
 !       p=w+alyam1*cO2i(k)
         p=w+alyam1*cO2i(k)+pgl(15,k,i,j)/pgl(4,k,i,j)
         h o2=const1*pgl(7,k,i,j)/g(k)
         x=(rads(k)+re)/h o2
         ch=chep(x,hi)
         pok=1.e-8*(pgl(1,k,i,j)*h o2*ch*1.e-4)**0.38
         if (pok.gt.60.) pok=60.
         p=p+alyam14*exp(-pok)
         pro= (rp(k)+rp(k-1))*.5
         a(k)=(1./dt+p)*pro
         b(k)=(q+cNo(k)/dt)*pro
       !  if(k.ne.nh) then
            dtn=pgl(7,k+1,i,j)-pgl(7,k-1,i,j)
            crab=0.5*dtn/pgl(7,k,i,j)/pro
        ! else
        !    dtn=pgl(7,k,i,j)-pgl(7,k-1,i,j)
        !    crab=dtn/pgl(7,k,i,j)/rp(k)
        ! end if
c     . . . scale height
         h =const*pgl(7,k,i,j)/g(k)
         sum=pgl(1,k,i,j)+pgl(2,k,i,j)+pgl(3,k,i,j)
         ams=(amo*pgl(3,k,i,j)+amo2*pgl(1,k,i,j)+
     *        amn2*pgl(2,k,i,j))/sum
         hsr=bk*pgl(7,k,i,j)/(ams*g(k))
c    . . . Coef. Mol. Dif.
         cmd(k)=3.e17/sum*sqrt(pgl(7,k,i,j))
         vrp=abs(pgl(10,k+1,i,j))
         vrm=abs(pgl(10,k-1,i,j))
         r(k)=cmd(k)*(1./h+crab)+ctd(k)*(1./hsr+crab)
  !       r(k)=r(k)-pgl(10,k,i,j)
         r(k)=r(k)-(pgl(10,k+1,i,j)+vrp)*.5-(pgl(10,k-1,i,j)-vrm)*.5
   1  continue
         a(1)=(rp(1)/dt)
         b(1)=(cNo(1)/dt)*rp(1)
         sum=pgl(1,1,i,j)+pgl(2,1,i,j)+pgl(3,1,i,j)
         cmd(1)=3.e17/sum*sqrt(pgl(7,1,i,j))
         dtn=pgl(7,2,i,j)-pgl(7,1,i,j)
         crab=dtn/pgl(7,1,i,j)/rp(1)
         h =const*pgl(7,1,i,j)/g(1)
         ams=(amo*pgl(3,1,i,j)+amo2*pgl(1,1,i,j)+
     *        amn2*pgl(2,1,i,j))/sum
         hsr=bk*pgl(7,1,i,j)/(ams*g(1))
         r(1)=cmd(1)*(1./h+crab)+ctd(1)*(1./hsr+crab)
         do k=1,nh-1
          cmd05=(cmd(k+1)+cmd(k)+ctd(k+1)+ctd(k))*.5
		r05=(r(k)+r(k+1))*.5
c          if(k.eq.nh-1) then
c            h =const*pgl(7,nh-1,i,j)/g(nh-1)
c            r05=cmd(nh-1)/h
c          end if
          d(k+1)=cmd05/rp(k)+r05*.5
          e(k+1)=(cmd05/rp(k)-r05*.5)/d(k+1)
          if(e(k+1).lt.0.) then
           d(k+1)=cmd05/rp(k)+r(k+1)
           e(k+1)=(cmd05/rp(k))/d(k+1)
          end if
        end do
        call potno(cNo,sdim,a,b,d,e,nh)
        do k=2,nh
         pgl(4,k,i,j)=cNo(k)
        end do
     	deallocate (a,b,d,e,r,cmd
     *          ,sdim)

	return
	end
 
!!

      subroutine baromi(an,pgl,rp,g,am,kpars,nh,its,ids,i,j,l)
      dimension an(nh)
     *         ,pgl(kpars,nh,its,ids),rp(nh),g(nh)
      data bk/1.38e-16/,ves/0.5/
      f2=bk/am
      do 3 k=l,nh
          tn=pgl(8,k-1,i,j)
          tv=pgl(8,k,i,j)
          h1=f2*tn/g(k-1)
          h2=f2*tv/g(k)
          alf=(ves/h1+(1.-ves)/h2)*rp(k-1)
          ss=alog(an(k-1)*tn/tv)-alf
          an(k)=exp(ss)
    3 continue
      return
	end
