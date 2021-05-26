!  O2 - 3-composition mixture                        23.12.2020
!  ver jan2020 2-D eddy diffusion coefficient   
      subroutine sdizkn_eddy(an1,an2,an3,an11,an21,an31,an6,vr,
     *                  vi,vj,ro,rp,r,g,n,n1,n2,dt,eddyco,ro1,
     *                  solu,gkoor,delta,nsu,dtet,uts,ddolg)
      USE mo_bas_gsm, ONLY:amo2,amn2,amo,qdis
      dimension an1(n1,n2,n),an2(n1,n2,n),an3(n1,n2,n)
     *         ,eddyco(n,n1),an11(n1,n2,n),an21(n1,n2,n),an31(n1,n2,n)
     *         ,vr(n1,n2,n),r(n),rp(n),g(n)
     *         ,ro(n1,n2,n),ro1(n1,n2,n)
     *         ,an6(n1,n2,n),vi(n1,n2,n),vj(n1,n2,n)
     *         ,solu(nsu),gkoor(2,n1,n2)
!     data amo2,amn2,amo/53.12e-24,46.51e-24,26.56e-24/
c
	allocatable q(:,:,:) ! dissosiation sourse
	allocate (q(n1,n2,n))
      n11=n1-1
      lk=17
!!!!!!!!!!!!!!!!!!!!!!! O2 altitude point !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      loov= 20   ! 20 - base
!     loov=28    ! for high activity

!!!!!!!!!!!!!!!!!!!!!!! O  altitude point !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     lov=21 ! winter var
!       lov=20  ! base var
       lov=26
c     lov=n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c    . . .  massive of fotodissociation rates (q)
      do i=1,n1
        do j=1,n2
          do k=1,n
		  q(i,j,k)=qdis(1,i,j,k)
	    end do
        end do
      end do

      n4=n
c         O    con,
      call gson_eddy(an1,an31,an2,an3,vr,an6,q,rp,r,g,lov,dt,
     *            eddyco,ro1,n1,n2,n ,vi,vj)
      call boskli(an31,n,n1,n2)
      call barsos(an31,an6,rp,g,amo,n,n1,n2,lov)
c . . . tridiagonal matrix algorithm 
!! comment 30.07.18 non explicit scheme
       call tn_jc(an31,vj,
     *            r,n,n1,n2,dt,n4)
       call boskli(an31,n,n1,n2)
       call tn_ic(an31,vi,
     *           r,n,n1,n2,dt,n4)
       call bongl(an31,n,n1,n2)
!!      call ficon(an31,r,n,n1,n2,dt,ddolg,vj,dtet,vi)
!!      call boskli(an31,n,n1,n2)
!!  end non explicit scheme 
!!          explicit scheme 
!      call tngojm(an31,vj,
!     *     r,n,n1,n2,dt,n)
!      call boskli(an31,n,n1,n2)
!      call tngoim_a(an31,vi,     ! advection across poles
!     *     r,n,n1,n2,dt,n)
!      call bongl(an31,n,n1,n2)
c . . . Old Var
c      call boskli(an31,n,n1,n2)
c      call barsos(an31,an6,rp,g,amo2,n,n1,n2,lov)
c . . . ï¢­ ï áå¥¬  ¯¥à¥­®á  ç¥à¥§ ¯®«îá
c       call tetcon(an31,r,n,n1,n2,dt,dtet,vi,vj,n4)
c      call tetcon_a(an31,vi,r,n,n1,n2,dt,n4)
c      call boskli(an31,n,n1,n2)
c . . .

c . . .   O2    con.
c
!       call gsoom (an1,an11,an2,an3,vr,an6,q,rp,r,g,loov,dt,
!     *              ctd,ro1,n1,n2,n ,vi,vj)
       call o2pro_3mix(an11,an1,an2,an3,an6,vr,vi,vj,
     *               q,eddyco,ro,r,rp,g,n1,n2,n,dt)
       call boskli(an11,n,n1,n2)

c  . . .   correction O2 above loov point
       call barsos(an11,an6,rp,g,amo2,n,n1,n2,loov)! loov)
!!   explicit scheme 
!       call tngojm(an11,vj,
!     *     r,n,n1,n2,dt,n)
!       call boskli(an11,n,n1,n2)
!       call tngoim_a(an11,vi,     ! advection across poles
!     *      r,n,n1,n2,dt,n) 
!       call bongl(an11,n,n1,n2)

!       call ficon(an11,r,n,n1,n2,dt,ddolg,vj,dtet,vi)
!       call boskli(an11,n,n1,n2)
c       call tetcon(an11,r,n,n1,n2,dt,dtet,vi,vj,n4)
c      call tetcon_a(an11,vi,r,n,n1,n2,dt,n4)
c      call boskli(an11,n,n1,n2)
c      call barsos(an11,an6,rp,g,amo2,n,n1,n2,loov)
c      call bongl(an11,n,n1,n2)
!! end of explicit scheme

c . . . tridiagonal matrix algorithm 
       call tn_jc(an11,vj,
     *            r,n,n1,n2,dt,n)  ! loov)
       call boskli(an11,n,n1,n2)
       call tn_ic(an11,vi,
     *            r,n,n1,n2,dt,n)  ! loov)
       call bongl(an11,n,n1,n2)
!       call boskli(an11,n,n1,n2)

c . . .   N2 con.
       key=0
c      lkk=lk
       lkk=n-1
       ot=amo2/amn2
       ot1=amo/amn2
       do 1 i=2,n11
        do 2 j=1,n2
         do 3 k=2,n-1
          tot=(ro1(i,j,k))/amn2
          sum=ot*an11(i,j,k)+ot1*an31(i,j,k)
          an21(i,j,k)=tot-sum
c         . . . check altitude point
          prov=0.5*an31(i,j,k)
c          prov=0.3*an31(i,j,k)
          if(an21(i,j,k).lt.prov) then
           if(k.lt.lkk) then
            key=1
            lkk= k-1
            if(lkk.LE.2) lkk=2 
            ii=i
            jj=j
           end if
          end if
    3    continue
    2   continue
    1  continue
       if (key.eq.1)
     *     print *,' N2 subtract until',lkk, ' point',ii,jj
c    *     print 100,an21(ii,jj,lkk),ii,jj,lkk
      call barsos   (an21,an6,rp,g,amn2,n,n1,n2,lkk)
c     call nts(an21,n,n1,n2,n4)
      call boskli(an21,n,n1,n2)
  100 format(' negative [N2]=',e10.2,3i4)
      
	deallocate (q)

      return
      end
