c     . . . Рассчет даты и времени
c
      subroutine dat  (ut,dayt,godt,verno,ut0,ut1,day,god)
      integer god,day,verno,dayt,godt
  900 format(' Ошибка !!! p/p dat:  (ut1-ut0) > utk   !!! '/
     *       ' ut1=',g10.2,'  ut0=',g10.2,'  utk=',g10.2)

      utk=3600*365*24
      if((ut1-ut0).le.utk)go to 4
        print 900,  ut1,ut0,utk
        verno=1
        go to 5
    4 continue
        dayt=day
        godt=god
        ut =ut1
        ng=god
        ngd=ng/4
        ngd=ng-ngd*4
        ngod=365
        if(ngd.eq.0)ngod=366
        i=0
    2   if(ut .lt.86400)go to 1
          ut =ut -86400
          i=i+1
          go to 2
    1   continue
        dayt=dayt+i
        if(dayt.le.ngod)go to 3
          dayt=dayt-ngod
          godt=godt+1
    3   continue
    5 continue
      return
      end
