      subroutine cpf4day(uts,day,nzap,ldor)
      allocatable pole(:)
      integer day
      character sr*2,srd*3
      character*80 fname
      allocate( pole(ldor/4))
      ihour=uts/3600.+0.001
	write(sr,'(i2)') ihour
        write(srd,'(i3)') day
	!print*,sr
	if(ihour.lt.10) sr(1:1)='0'
      if(day.lt.100) srd(1:1)='0'
      if(day.lt.10) srd(1:2)='00'
      fname='file4_'//srd//'.'//sr  
      open(44,file=(fname),
     *       access='direct',form='unformatted',recl=ldor) 
        do irecl=1,nzap
   !        print*,'irecl=',irecl
           read(4,rec=irecl) pole
           write(44,rec=irecl) pole
        end do
        close(44)
        deallocate(pole)
	return
	end