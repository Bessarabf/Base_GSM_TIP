c cyclt2-vdrift-bdip ;
c        trubka-tube-bdip;
c                    lambet-bdip;
c                    vorvfv-bdip;
c                    backfv-bdip;
c                    ntnb-bdip;
c                    alga-bdip;
c
!      function bdip(alt,tet)
!      double precision d,t
!      data re/6371.02e5/,g10/-.30356/
!      t=tet
!      ro=1./(1.+alt/re)
!  900 format(' ',10g12.4)
!      roq=ro**3
!      d=dsqrt(1.d0+3.d0*dcos(t)**2)
!      bdip=-g10*roq*d
!      return
!      end
! Updated Mihail Melnik from Polar Geophysical Institute 14.02.2020
!        function bdip(alt,tet)
!         real,parameter :: re = 6371.02e5, g10 = .30356
!         real ro,alt,tet
!         real*8 t
!         ro = re/(re+alt)
!         t = dcos(dble(tet)) * dcos(dble(tet))
!         bdip = g10 * ro * ro * ro * real(dsqrt(1.d0+3.d0*t ))
!        end
! Updated Mihail Melnik from Polar Geophysical Institute 06.05.2021
       function bdip(alt,tet)
        real,parameter :: re = 6371.02e5, g10 = .30356
        real ro,alt,tet
        real t
        ro = re / ( re + alt )
        t  = cos( tet ) * cos( tet )
        bdip = g10 * ro * ro * ro * sqrt( 1.0 + 3.0 * t )
       end
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