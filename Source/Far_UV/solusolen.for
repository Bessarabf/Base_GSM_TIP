!!!!!!!  solusolen uses in var GSM TIP with modeldan 
      subroutine solusolen(solu,sole,solen,solet,nsu,nse,mas)     
      dimension solu(nsu),sole(nse),solen(nse),solet(nse),mas(*)
 !    
      integer god
      character*60,fname,smas*2
!!! form file''s name 
      Write(smas, '(i2)' )  mas(4)
      fname='solen'//smas//'.dat'
      !open(22,file=trim(fname),status='old')
	open(22,file=(fname),status='old')

      read(22,*)
      read(22,*) fa0,fs
      nn=mas(4)
      print*,'mas(8)=',mas(8),nse
      if(mas(8).ne.0) then
          if (mas(8).eq.1 ) then
     	     call flosu(god,fa0,sole,nn)
	     print*,' flosu'
	    else if (mas(8).eq.2) then
              if(mas(4).le.15) then
     	         call flosuN(fs,fa0,sole,nn)         ! Default option 
	         print*,' flosuN'
              else
                 call flosuN38(fs,fa0,sole,mas(4))   ! Option for 38 bands (and flares)
                 print*,' flosuN38'
	        end if
          else if (mas(8).eq.3) then
              if(mas(4).le.15) then
                 call flosuEUVAC(fs,fa0,sole,nn)
                 print*,' flosuEUVAC'
              else
                 call flosuEUVAC39(fs,fa0,sole,mas(4))
                 print*,' flosuEUVAC39'
	        end if
          else
            print*,' wwod: incorrect mas(8)'
            stop
 	    end if
          
!!! Nusinov solu model (default option)
          !!!! La in solu 
          solu(1)=sole(nse) ! La in photon
          solu(nsu/2+1)=solu(1)*1.98648/121.5 !La in erg
          call nus_uv(solu,nsu)
      end if
      read(22,*)
!!!! reading solen & sole (if mas(8) = 0)      
      do i=1,nse
         if (mas(8).eq.0) then
             read(22,*,err=100) dws,dwe, sole(i),solet(i),solen(i)
!!!!             print*,i
         else 
             read(22,*,err=100) dws,dwe, tmp,solet(i),solen(i)
         end if
      end do
!!!! read solu (if nessesary)
      if (mas(8).eq.0) then    
	read(22,*)
	read(22,*)
        do i=1,nsu/2     
           read(22,*,err=100) dws,dwe,solu(i),solu(i+nsu/2)
        end do
!!!! correct La in solu 
        solu(1)=sole(nse) ! La in photon
        solu(nsu/2+1)=solu(1)*1.98648/121.5 !La in erg
      end if 
      close(2)
      return
  100 print*, ' incorrect nsu or nse',nsu,nse
      stop     
      end	
