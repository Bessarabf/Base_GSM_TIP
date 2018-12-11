c                Расчет даты и времени
      subroutine calcdat(ut0,ut1,day,god,ng,nden,namemes,np,min)
      integer god,day
      integer kdm(12)/31,29,31,30,31,30,31,31,30,31,30,31/
      character*8 mes(12),namemes
      data mes/'january ','february','march   ','april   ',
     *         'may     ','june    ','july    ','august  ',
     *         'septembr',
     *         'october ','november','december'/
  920 format(i19,' year',i3,1x,a8,i3,' hour.',i2,' min.'/)
c                 Расчет времени суток, дня и года
      nd=day
      ut=ut1
      ng=god
      ngd=ng/4
      ngd=ng-ngd*4
      ngod=365
      if(ngd.eq.0)ngod=366
      i=0
    1 if(ut.lt.86400)go to 2
        ut=ut-86400
        i=i+1
        go to 1
    2 continue
      nd=nd+i
      if(nd.le.ngod)go to 3
        nd=nd-ngod
        ng=ng+1
    3 continue
c                 Расчет дня месяца и номера месяца
      np=ng/4
      np1=ng-np*4
      if(np1.eq.0)go to 4
        kdm(2)=28
    4 continue
      nmes=1
      nden=nd
      if(nden.le.31)go to 7
        do 6 i=1,12
          nden=nden-kdm(i)
          nmes=nmes+1
          if(nden.gt.kdm(i+1))go to 5
            go to 7
    5     continue
    6   continue
    7 continue
c          nmes - номер месяца
c          nden - день месяца
      nsek=0
      np=0
      min=0
      nut0=ut
      np=nut0/3600
      np1=nut0-np*3600
      min=np1/60
c     nsek=np1-min*60
      namemes=mes(nmes)
c     print   920,ng,nden,namemes,np,min
      return
      end
