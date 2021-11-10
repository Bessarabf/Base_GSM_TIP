      function oml(dolm)
      data st0/1.9765732e-1/,om/7.272205e-5/,
     *pi/3.14159265359/
      fi=dolm/180.*pi
      oml=om*st0*sin(fi)
      return
      end
