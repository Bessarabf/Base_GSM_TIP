! Updated Mihail Melnik from Polar Geophysical Institute 06.05.2021
      function pite(cio,cih,cihe,te)
!      data c1/3.70e-3/,c2/5.92e-2/,c3/1.48e-2/,s/-1.5/
!      real,parameter :: c1 = 3.70e-3, c2 = 5.92e-2, c3 = 1.48e-2, s = -1.5
      real,parameter :: c1 = 3.70e-3, c2 = 5.92e-2, c3 = 1.48e-2
!     pite=(c1*cio+c2*cih+c3*cihe)*te**s
      pite= ( c1 * cio + c2 * cih + c3 * cihe ) / sqrt( te * te * te )
      return
      end
