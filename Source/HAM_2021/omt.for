c   omt :  gst, lambet, 
c
      function omt(tet,dolm)
      data ct0/9.8027118e-1/,st0/1.9765732e-1/,om/7.272205e-5/,
     *pi/3.14159265359/
      fi=dolm/180.*pi
      omt=-om*(sin(tet)*ct0+cos(tet)*st0*cos(fi))
      return
      end
