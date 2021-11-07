      subroutine inter8(plm1,plm2,plm3,kpart,uo,du,
     *                  dolo,ddolgt,ux,dolx,log,nv,l,pm,plm8)
	  
      dimension plm1(kpart),plm2(kpart),plm3(kpart),del(3)
      dimension pm(kpart,nv),plm8(kpart)
      data k/1/,al/0.54e-10/
      
      IF(ABS(DU).LT.(1.E-6))THEN
	      PRINT *,'INTER8 ',UO,DU,L,DOLO,DDOLGT,DOLX
!		 PAUSE
		 DU=1.E-6
		 dux=0.
      else
        dux=(ux-uo)/du
      END IF



      ddolx=(dolx-dolo)/ddolgt
      do 1 i=1,KPART !    15.04.15 ! I=1,8
        p1=plm1(i)
        p2=plm2(i)
        p3=plm3(i)
        p8=plm8(i)
        if(i.gt.3)go to 2
          if(p1.lt.al)p1=al
          if(p2.lt.al)p2=al
          if(p3.lt.al)p3=al
          if(p8.lt.al)p8=al
          if(log.eq.0)go to 3
            p1=alog(p1)
            p2=alog(p2)
            p3=alog(p3)
            p8=alog(p8)
   3      continue
   2    continue
        p4=p1+(p2-p1)*dux
        p5=p3+(p8-p3)*dux
        pm(i,l)=p4+(p5-p4)*ddolx

    1 continue
      return
      end
