c    . . . Вариант начала счета (определяется UPRS & UPRT)
      subroutine nachus(uprt,uprs,ut0,nh,dtets,dtett,ddolgs,ddolgt,
     *                 day,ap0,fa0,fs,ntsl,nl,god,ldor,nzapt,nzaps,
     *                 nadrt,nadrs,kdf,kpars,kpart,rads,park,ks,its,
     *                 kdu,isp,par,pole,nr,verno,pkp0,mass,imja,mast)
      dimension ntsl(nl),par(nr),pole(ldor/4),rads(nh),mast(40),
     *          kdf(20),kdu(20),park(ks),mass(30)
      character *1 imja(80)
      integer uprt,uprs,day,god,verno
      logical readfl
      mdor=ldor/4
c                         Shar
      nfw=11
      nfr=5
      nfile=5
      if(uprs.eq.0)then ! нулевые начальные условмя - не применяются
c
c        call zeros(day,ap0,fa0,fs,ut0,nh,dtets,ddolgs,rads,
c     *       ldor,kdf,isp,kpars,par,pole,nr,verno,its,pkp0,mass)
        print   921
  921   format(' zeros(shar) -> Nachus  ')
c       write (10,921)
        readfl=.false.
        call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *            god,day,ut0,pole,verno,imja)
      else if(uprs.eq.2) then ! время счета задается в danmodel 
c           write (10,941)
            print   941
  941       format(' Disk(shar) -> Nachus   ')
            readfl=.true.  ! читаем текущее время и дату в file4
            call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *                nt,nt1,unt,pole,verno,imja)
            readfl=.false. ! пишем дату и время   
            call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *                god,day,ut0,pole,verno,imja)
      else if(uprs.eq.3) then ! продолжение счета
    
            print*, ' Calculation;   uprt = uprs = 3 '
            readfl=.true.  ! читаем дату и время 
            call rtime(readfl,nfile,nzaps,nadrs,isp,ldor,kdf,uprt,
     *                nt,nt1,unt,pole,verno,imja)
           ! go to 10
      else 
            print*, ' Incorrect uprs -> stop !   '
            stop      
      end if
      if(verno.eq.1) return
c
c                       Trubka
      nfw=12
      nfr=6
      nfile=6
      readfl=.false.
      if(uprt.eq.0) then
            call zerot(ddolgt,ntsl,nl,ut0,kdf,ldor,isp,kpart,
     *      par,pole,nr,park,ks,mast)
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *       god,day,ut0,pole,verno,imja)
            print*,' zerot(trubka) -> Nachus '
      else if(uprt.eq.1) then 
c           call mldisk(nfile,kdf,kdu,ldor,isp,par,nr,verno)
            readfl=.true.
            print*, ' ML(trubka) ->  Nachus '
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 nt,nt1,unt,pole,verno,imja)
            readfl=.false.
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 god,day,ut0,pole,verno,imja)
      else if(uprt.eq.2) then 
            print*, ' Disk(trubka) -> Nachus '
            readfl=.true.
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 nt,nt1,unt,pole,verno,imja)
            readfl=.false.
            call rtime(readfl,nfile,nzapt,nadrt,isp,ldor,kdf,uprt,
     *                 god,day,ut0,pole,verno,imja)
      end if
      return
      end
