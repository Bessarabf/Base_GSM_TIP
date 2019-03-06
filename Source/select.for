!     VER 18.04.14 add nv to interface      
	subroutine select(ne,j,m,ntsl,nl,cio1,cih1,cihe1,vio1,vih1,
     *           vihe1,ti1,te1,co2,cn2,co,ch,che,cim,tn,vnq,vnu,vnv,
     *           qo,qsm,ht,tt,vdv,vdu,par,nr,NV)

      dimension ntsl(nl),par(nr),cio1(nv),cih1(nv),cihe1(nv),vio1(nv),
     *          qsm(nv),
     *          vih1(nv),vihe1(nv),ti1(nv),te1(nv),co2(nv),cn2(nv),
     *          co(nv),ch(nv),che(nv),cim(nv),tn(nv),
     *          vnq(nv),vnu(nv),vnv(nv),qo(nv),ht(nv),tt(nv),
     *          vdv(nv),vdu(nv)
      data pi/3.14159265/
  900 format(' ',10g12.4)
      cr=pi/180.
      i1=0
      if(j.NE.1) THEN
        i2=j-1
        do i=1,i2
          i1=i1+ntsl(i)
        end do
      END IF
      n=ntsl(j)
      do 7 i=1,n
        k=(i1+i-1)*m+1
        goto(3,4,5,6),ne
    3   continue
          cio1(i)=par(k)
          cih1(i)=par(k+1)
          cihe1(i)=par(k+2)
          vio1(i)=par(k+3)
          vih1(i)=par(k+4)
          vihe1(i)=par(k+5)
          ti1(i)=par(k+6)
          te1(i)=par(k+7)
!	print*,n,i,ti1(i),te1(i)
          goto7
    4   continue
          co2(i)=par(k)
          cn2(i)=par(k+1)
          co(i)=par(k+2)
          cim(i)=par(k+3)
          qo(i)=par(k+4)
          ch(i)=par(k+5)
          che(i)=par(k+6)
          tn(i)=par(k+7)
          vnq(i)=par(k+8)
          vnu(i)=par(k+9)
          vnv(i)=par(k+10)
          qsm(i)=par(k+11)
          
          goto7
    5   continue
          ht(i)=par(k)
          tt(i)=par(k+1)*cr
          goto7
    6   continue
          vdv(i)=par(k)
          vdu(i)=par(k+1)
    7 continue
      return
      end
