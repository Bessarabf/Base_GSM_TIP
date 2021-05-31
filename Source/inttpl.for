!      ver. 20.04.14 
	subroutine inttpl(mi,dolg,dolg1,vv,par,par1,nr,nfile,kpart,
     *           ddolgt,ntsl,nl,kdf,POLE,ldor,isp,mast,dolga,dolgb)
      dimension par(nr),par1(nr),ntsl(nl),kdf(20)
	
      dimension mast(40),pole(ldor/4)
      logical readfl
      	
      readfl=.true.
        k1=3
        k2=3
        md=1
        if(mi.ne.1) go to 1
          call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *    kdf,ldor,isp,md,par,pole,nr,mast)
          dolga=dolg
          if(vv.eq.0.) go to 2
            call wwt(readfl,nfile,kpart,dolg1,ddolgt,1,ntsl,nl,
     *      kdf,ldor,isp,md,par1,pole,nr,mast)
            dolgb=dolg1
            go to 3
    2     continue
            do 4 i=1,nr
              par1(i)=par(i)
    4       continue
            dolgb=dolg
    3     continue
          mi=2
          go to 22
    1   continue
          if(vv.eq.0.) go to 8
            if(dolg1.eq.dolgb) go to 7
              k2=0
              if(dolg1.eq.dolga) k2=1
!              print *,'dolgb=',dolgb,'dolg1=',dolg1
              dolgb=dolg1
    7       continue
            go to 9
    8     continue
            k2=2
    9     continue
          IF(dolg.NE.dolga) THEN
            k1=0
            if(dolg.eq.dolgb) k1=1
            dolga=dolg
          END IF
          if(vv.eq.0.) dolgb=dolg
          IF(k1.EQ.0) THEN
            if(k2.EQ.0) THEN
              call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par,pole,nr,mast)
              call wwt(readfl,nfile,kpart,dolg1,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par1,pole,nr,mast)
            end if
            IF(k2.EQ.1) then
              !!!!!!!!!!!!!!!!!!!do 13 i=1,nr
                par1=par
              !!!!!!!!!!!!!!!!!!!13         continue
              call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par,pole,nr,mast)
            END IF
            if(k2.EQ.2) THEN
              call wwt(readfl,nfile,kpart,dolg,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par,pole,nr,mast)
              !!!!!!!!!!!!!!!! do 14 i=1,nr
                par1= par !!!!!!!!!!par1(i)=par(i)
              !!!!!!!!!!!!!!!!!!14         continue
            END IF
           END IF
          if(k1.EQ.1) then
            if(k2.EQ.0) then
            !!!!!!!!!!!!!!!!!!!!!!          do 17 i=1,nr
              par=par1                     !!!par(i)=par1(i)
            !!!!!!!!!!!!!!!!!!!!!!!!!! 17         continue
              call wwt(readfl,nfile,kpart,dolg1,ddolgt,1,ntsl,nl,
     *        kdf,ldor,isp,md,par1,pole,nr,mast)
            end if  
            if(k2.eq.1) print 900
  900       format(' inttpl:    incorrect k2 ')
            if(k2.EQ.2) then
            !!!!!!!!!!!! do i=1,nr
                par=par1 !!!!!!!!!par(i)=par1(i)
            !!!!!!!!!!!!!19         continue
            end if
          end if
   22   continue
cc      mi=2


        return
        end

