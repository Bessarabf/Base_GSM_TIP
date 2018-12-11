      subroutine vmion(pgl1,kpars,nh,its,rads,dtett,
     *           dtet,ddolgs,potef,ntr,nl2,ids,vim,vid,vir,mast,mass)
      dimension pgl1(kpars,nh,its,ids),mast(40),mass(30),rads(nh),
     *          vim(nh,its,ids),vid(nh,its,ids),vir(nh,its,ids),
     *          potef(ntr,ids,nl2)
      real nu0
      double precision ami,ae,ame,cr
      data nu0/1.e-9/,ami/30./,ae/1.6e-24/,e/1.6e-20/,b0/0.3/,
     *     pi/3.1415926/,re/6371.02e5/
      cr=180./pi
      ame=ae*ami/e
      nmd=360./ddolgs+0.1
      iteq=(its+1)/2
      itd=its-1
      do 12 j = 1 , ids
        do 11 itet = 2 , itd
c. . . двойной шаг для 10 - 5 гр. сетки
c        it=itet*2-1
c
      tet=dtet*(itet-1)/cr
      dip=1.
      if(tet.eq.0.)dip=pi/2.
      if(tet.eq.pi)dip=-pi/2.
      if(tet.eq.pi/2.)dip=0.
      if(dip.ne.1.)go to 5
        t=tan(tet)
        dip=atan(2./t)
    5 continue
      si=sin(dip)
      ci=cos(dip)
      sc=si*ci
      ct=cos(tet)
      st=sin(tet)
      bb=b0*sqrt(1+3.*ct*ct)
      be=ame/bb
      bnu=be*nu0
c
c 774 format (' vmion - 2 ')
c     print 774
      dd2=2*ddolgs/cr
      if(j.ne.1)go to 2
        def=(potef(ntr,2,itet)-potef(ntr,nmd,itet))/dd2
        go to 3
    2 continue
      if(j.ne.nmd)go to 4
         def=(potef(ntr,1,itet)-potef(ntr,nmd-1,itet))/dd2
         go to 3
    4 continue
c 775 format (' vmion - 3 ')
c     print 775
      def=(potef(ntr,j+1,itet)-potef(ntr,j-1,itet))/dd2
    3 continue
      defs=def/st/bb
      det=(potef(ntr,j,itet+1)-potef(ntr,j,itet-1))/2./dtett*cr
c . . . Граница замкнутых-разомкнутых
      if(itet.eq.4.or.itet.eq.33)det=(potef(ntr,j,itet)
     *                           -potef(ntr,j,itet-1))/dtett*cr
      if(itet.eq.5.or.itet.eq.34)det=(potef(ntr,j,itet+1)
     *                           -potef(ntr,j,itet))/dtett*cr
      detb=det/bb
      r=re+rads(ntr)
      et=-detb/r
      el=-defs/r
c
      if(itet.ne.iteq) go to 6
        tetd=tet-dtett/cr
        stet=1.-sin(tetd)**2
        er=-(potef(ntr,j,itet)-potef(ntr,j,itet-1))/(r*stet)
        er=er/bb
        go to 7
    6 continue
      er=-et*ci/si
    7 continue
      do 1 kk=1,nh
        g=bnu*((pgl1(1,kk,itet,j)+pgl1(2,kk,itet,j))/2.
     * +pgl1(3,kk,itet,j)/3.)
        g2=g*g
        alf=1.+g2
        etb=et/g
        etl=el/g
        etr=er/g
        a=pgl1(10,kk,itet,j)+etr
        b=pgl1(11,kk,itet,j)+etb
        c=pgl1(12,kk,itet,j)+etl
c       *****
        if (mast(13).eq.0) then
        a=pgl1(10,kk,itet,j)
        b=pgl1(11,kk,itet,j)
        c=pgl1(12,kk,itet,j)
        end if
c       *******
        if (mass(6).eq.0) then
        a=etr
        b=etb
        c=etl
        end if
c       *******
c 678 format (' vmion - cycl')
c     print 678
        vir(kk,itet,j)=((si*si+g2)*a+sc*b+g*ci*c)/alf
        vim(kk,itet,j)=((ci*ci+g2)*b+sc*a-g*si*c)/alf
        vid(kk,itet,j)=(g*si*b-g*ci*a+g2*c)/alf

    1 continue
   11 continue
   12 continue
	! ╙тхышўхэшх ьюы. ёъюЁюёЄш т Єюўърї:
!	do kk=1,nh
! 	      vir(kk,8,2)=10*vir(kk,8,2)
!		  vir(kk,9,2)=10*vir(kk,9,2)
!	      vir(kk,8,3)=10*vir(kk,8,3)!
!	      vir(kk,9,3)=10*vir(kk,9,3)
!	
!         	vim(kk,8,2)=10*vim(kk,8,2)
!		vim(kk,9,2)=10*vim(kk,9,2)
!	    vim(kk,8,3)=10*vim(kk,8,3)
!	    vim(kk,9,3)=10*vim(kk,9,3)
!
!         	vid(kk,8,2)=10*vid(kk,8,2)
!		vid(kk,9,2)=10*vid(kk,9,2)
!	    vid(kk,8,3)=10*vid(kk,8,3)
!	    vid(kk,9,3)=10*vid(kk,9,3)
!	end do
c 900 format(' ',10g12.3)
c 901   format(' vir ')
c 902   format(' vim ')
c 903   format(' vid ')
c 777 format (' vmion - wihod ')
c     print 777
      return
      end
