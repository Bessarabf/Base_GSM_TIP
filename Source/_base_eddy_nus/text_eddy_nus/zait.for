c . . . Запись и чтение времени и даты
      subroutine zait(readfl,nfile,nzapt,nzaps,nadrt,nadrs,isp,
     *                ldor,kdf,god,day,ut0,ut1,pole)
      integer nfa(20)/4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,
     *       0,0,0,0,0/
      dimension kdf(20),pole(ldor/4)
      logical readfl
      integer adres,god,day
      mdor=ldor/4
      nzap=nzapt
      adres=nadrt
      if(nfile.eq.6)go to 4
        nzap=nzaps
        adres=nadrs
    4 continue
      isp=kdf(nfile)+nzap
      nf=nfa(nfile)
      read(nf,rec=isp)(pole(i),i=1,mdor)
      if(readfl) go to 1
        pole(adres+3)=god
        pole(adres+2)=day
        pole(adres+1)=ut1
        pole(adres)=ut0
        isp=kdf(nfile)+nzap
        write(nf,rec=isp)(pole(i),i=1,mdor)
        go to 2
    1 continue
        ut1=pole(adres+1)
    2 continue
c      print 900,ut0,ut1,god,day,readfl
      return
      end

