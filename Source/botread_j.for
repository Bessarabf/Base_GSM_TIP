c         Считывание данных Jacobi etc.
      subroutine botread(day,pril,its,ids,kpa,nt) ! add parameter day
                                                  ! day - month day
      dimension pril(kpa,its,ids,nt)
      integer day
      character*2 cday,str*72

      write(cday,'(I2)') day

	print*,'botread cday=',cday
	
      if(day.lt.10)cday(1:1)='0'

      open(9,file='HAMMONIA/PL_'//cday,status='old')
      open(10,file='HAMMONIA/Vwind_'//cday,status='old')
      open(11,file='HAMMONIA/Zwind_'//cday,status='old')
      open(12,file='HAMMONIA/T_'//cday,status='old')
c
     
      do l = 1,nt
       DO i = 1 , ITS
        read(9,*)str

         READ (9,*) (pril(1,i,j,l),j=1,ids)
     !
         read(12,*)str
         READ (12,*) (pril(2,i,j,l),j=1,ids) !Tn
      !
         read(10,*)str
         READ (10,*) (pril(4,i,j,l),j=1,ids) ! meridional
      !
         read(11,*)str
         READ (11,*) (pril(5,i,j,l),j=1,ids)
      enddo
	

      enddo
      do l = 1,nt
       do j=1,ids
        do i=1,its
            pril(3,i,j,l)=0.
            pril(4,i,j,l)=-100.*pril(4,i,j,l)
            pril(5,i,j,l)=100.*pril(5,i,j,l)
        enddo
      enddo

      enddo

  900 format (a25)
  901 format (a1)
  902 format (F5.1,4F8.2)
  903 format (F5.1,3X,E9.2,8X,E9.2,5X,3F8.1)
  904 format (5E11.3)
  101 format(1x,2x,6e12.4)
      
      close(9)
      close(10)
      close(11)
      close(12)
      return
      end

