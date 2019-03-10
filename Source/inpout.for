!
! чтение - запись параметров трубки
! readfl - условие чтения - записи
! kdf - массив, определяющий начальную запись для чтения
! ldor - длина записи (байт)
! isp - номер записи
! nob - 1
! lpar - 
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
      if(nfile.gt.8)nf=5
      ndor=nob/mdor
      nost=nob-ndor*mdor
      nsvob=mdor-nost
      if(nost.eq.0)nsvob=0
      if(lpar.gt.nsvob)go to 18
        kd=0
        ln=0
        kz=0
        lost=0
        go to 17
   18 continue
      ln=lpar-nsvob
      kd=ln/mdor
      kz=kd*mdor
      lost=ln-kz
   17 continue
      isp=ndor+kdf(nfile)+1
      if(readfl) go to 7 ! идти читать, иначе запись
        if(nsvob.eq.0)go to 2
          read(nf,rec=isp)(pole(i),i=1,mdor)
          n=nost+1
          k=n+lpar-1
          if(nsvob.le.lpar)k=mdor
          i=1
          do 1 j=n,k
            pole(j)=par(i)
            i=i+1
    1     continue
          write(nf,rec=isp)(pole(i),i=1,mdor)
    2   continue
        if(kd.eq.0) go to 3
          n=nsvob+1
          k=n+mdor-1
          is=isp
          if(nsvob.ne.0)is=is+1
          do 22 m=1,kd
            write(nf,rec=is)(par(i),i=n,k)
            n=n+mdor
            k=k+mdor
            is=is+1
   22     continue
    3   continue
        if(lost.eq.0)go to 5
          is=isp
          if(nsvob.ne.0)is=isp+1
          if(kd.ne.0)is=is+kd
!          print*,'inpouT',' rec=',isp
          read(nf,rec=is)(pole(i),i=1,mdor)
          i=1
          n=nsvob+kz+1
          k=n+lost-1
          do 4 j=n,k
            pole(i)=par(j)
            i=i+1
    4     continue
          write(nf,rec=is)(pole(i),i=1,mdor)
    5   continue
        go to 6
    7 continue ! чтение
        if(nsvob.eq.0)go to 8
          read(nf,rec=isp)(pole(i),i=1,mdor)
          isp=isp+1
          k=lpar
          if( nsvob.le.lpar)k=nsvob
          i=nost+1
          do 10 j=1,k
          par(j)=pole(i)
            i=i+1
   10     continue
    8  continue
        if(ln.eq.0)go to 12
          n=nsvob+1
          k=ln+nsvob
          k1=n+mdor-1
          if(ln.le.mdor) k1=k
          m=kd
          if(lost.ne.0)m=m+1
          do 11 j=1,m
            read( nf,rec=isp)(par(i),i=n,k1)
            isp=isp+1
            n=n+mdor
            k1=k1+mdor
            if(k1.gt.k)k1=k
  11      continue
   12   continue
    6 continue
      return
      end

