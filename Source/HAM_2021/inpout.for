!
! чтение - запись параметров трубки
! readfl - условие чтения - записи
! kdf - массив, определяющий начальную запись для чтения
! ldor - длина записи (байт)
! isp - номер записи
! nob - 1
! nf  - номер файла в open(nf
! nfile - номер элемента в kdf - начало записи shar, pot и т.д.
! lpar = nh*kpar 
! par - одномерный массив с параметрами, который заполняется
! pole - рабочий массив, сюда читаем данные
! nr -  
      subroutine inpout(readfl,nfile,kdf,ldor,
     *       isp,nob,lpar,par,pole,nr)
      logical readfl
      dimension par(nr),pole(ldor/4)
      dimension kdf(20)
!	if (nfile.eq.5) then
!	 print*,kdf
!	 stop
!	end if
      mdor=ldor/4
      nf=4
      if(nfile.gt.8) then 
         nf=5 ! файл f5
      !print*, ' inpout file5', readfl
	end if
!    integer number of records
      ndor=nob/mdor
!    difference of addresses
      nost=nob-ndor*mdor
!    
      nsvob=mdor-nost
      if(nost.eq.0)nsvob=0
      if(lpar.le.nsvob) then
        kd=0
        ln=0
        kz=0
        lost=0
      else
        ln=lpar-nsvob
        kd=ln/mdor
        kz=kd*mdor
        lost=ln-kz
      end if
      isp=ndor+kdf(nfile)+1
      if(readfl) go to 7 ! идти читать, иначе запись
! 
        if(nsvob.ne.0) then 
          read(nf,rec=isp) pole
          n=nost+1
          k=n+lpar-1
          if(nsvob.le.lpar)k=mdor
          i=1
          do j=n,k
            pole(j)=par(i)
            i=i+1
          end do 
          write(nf,rec=isp) pole 
        end if 
        if(kd.ne.0) then 
          n=nsvob+1
          k=n+mdor-1
          is=isp
          if(nsvob.ne.0)is=is+1
          do m=1,kd
            write(nf,rec=is)(par(i),i=n,k)
            n=n+mdor
            k=k+mdor
            is=is+1
          end do 
        end if 
        if(lost.ne.0) then
          is=isp
          if(nsvob.ne.0)is=isp+1
          if(kd.ne.0)is=is+kd
!          print*,'inpouT',' rec=',isp
          read(nf,rec=is) pole
          i=1
          n=nsvob+kz+1
          k=n+lost-1
          do j=n,k
            pole(i)=par(j)
            i=i+1
          end do
          write(nf,rec=is) pole
        end if 
        go to 6
    7 continue ! чтение
        if(nsvob.ne.0) then
          read(nf,rec=isp) pole
          isp=isp+1
          k=lpar
          if( nsvob.le.lpar) k=nsvob
          i=nost+1
          do j=1,k
            par(j)=pole(i)
            i=i+1
          end do
        end if 
        if(ln.ne.0) then
          n=nsvob+1
          k=ln+nsvob
          k1=n+mdor-1
          if(ln.le.mdor) k1=k
          m=kd
          if(lost.ne.0)m=m+1
          do j=1,m
            read( nf,rec=isp)(par(i),i=n,k1)
            isp=isp+1
            n=n+mdor
            k1=k1+mdor
            if(k1.gt.k)k1=k
          end do 
        end if
    6 continue
      return
      end

