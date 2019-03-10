c              Чтение или запись даты и времени
c     ver 29/05/2014
 
      subroutine rtime(readfl,nfile,nzap,adres,isp,ldor,kdf,uprt,
     *                god,day,ut0,pole,verno,imja)
      logical readfl

      integer uprt,verno,god,day,adres

      character*1 ntexto(80),imja(80)
      dimension pole(ldor/4),imjan(80),imjar(80),kdf(20)

      integer nfa(20)/4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,
     *       0,0,0,0,0/

      integer kdm(12)/31,28,31,30,31,30,31,31,30,31,30,31/
      character*8 mes(12)
      data mes/'january ','february','march   ','april   ',
     *         'may     ','june    ','july    ','august  ',
     *         'septembr','october ','november','december'/
      mdor=ldor/4

      isp=kdf(nfile)+nzap

      nf=nfa(nfile)
	print*,nf,isp,nfile,kdf(nfile),nzap
	
      read(nf,rec=isp)pole
 	print*,nf,isp,pole(adres+3)
	
      if(.not.readfl) then ! write data 
         pole(adres+3)=god
         pole(adres+2)=day
         pole(adres)=ut0
         pole(adres+1)=ut0

! . . . преобразуем название месяца 
         call kodir(imja,imjan,0,imjar,ntexto)
         jjj=3
         do iii=1,80
            pole(adres+jjj+iii)=imjan(iii)
         end do
! пишем название месяца за датой
        isp=kdf(nfile)+nzap
        write(nf,rec=isp) pole
      else  ! читаем дату и время  go to 8
    
        ut0=pole(adres)
        ut1=pole(adres+1)
        day=pole(adres+2)
        god=pole(adres+3)

        jjj=3

        do iii=1,80
          imjar(iii)=pole(jjj+iii+adres)
        end do 
! . . .преобразуем в текст
        call kodir(imja,imjan,1,imjar,ntexto)
        if(uprt.ge.2)then
          ut=ut0
          nd=day
          ng=god
          call dat(ut0,day,god,verno,ut,ut1,nd,ng)
        end if
     
        if(verno.eq.1) return ! end of rtime
! . . . проверка на високосность          
          np=god/4
          np1=god-np*4
          kdm(2)=29
          if(np1.ne.0) kdm(2)=28 ! february
     
          nmes=1
          nden=day
! . . . ищем месяц и день месяца:
          if(nden.gt.31) then 
          
			do while(nden.gt.kdm(nmes)) 
				nden=nden-kdm(nmes)
				nmes=nmes+1
				if(nmes.gt.12) then 
				  print*, ' incorrect number day=',day,' STOP'
				  stop
				end if 
			end do 
	    end if
c     nmes- nomer mesjaca
c     nden- den  mesjaca
          nsek=0
          np=0
          min=0
          nut0=ut0
          np=nut0/3600
          np1=nut0-np*3600
          min=np1/60
          nsek=np1-min*60
          print   920,god,nden,mes(nmes),np,min,nsek,ntexto
    9   continue
      end if

  920 format(' *',i5,' year',i3,1x,a8,' time: ',
     *   i3,' hour.',i2,' min.',i2,' sek. ,',8a1/' ',72a1)

      return
      end

