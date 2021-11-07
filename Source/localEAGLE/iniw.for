      subroutine iniw(ns,wm)
      dimension wm(*)
      dimension gis(24)
      data gis/1.,60.,120.,165.,180.,205.,280.,310.,460.,630.,800.,
     *         911.,1037.,1191.,1386.,1632.,1955.,2389.,2998.,3897.,
     *	   5326.,7874.,13459.,33703./
      do i=2,ns+1
        w1=12397./gis(i-1)
        w2=12397./gis(i)
	  wm(i-1)=(w1+w2)/2.
      end do
      return
      end