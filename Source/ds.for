!      function ds(altm,alt,altp,tetm,tetp)
!      data re/6371.02e5/
!      ds=sqrt((altp-altm)**2+(alt+re)**2*(tetp-tetm)**2)
!      return
!      end
! Updated Mihail Melnik from Polar Geophysical Institute 14.02.2020
       function ds(altm,alt,altp,tetm,tetp)
       real,parameter :: re = 6371.02e5
        real altm,alt,altp,tetm,tetp
       ds=sqrt((altp-altm)*(altp-altm)+(alt+re)*(alt+re)*(tetp-tetm)*
     *  (tetp-tetm))
       end