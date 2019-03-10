      function del(day)
      integer day
      data tg/4.3481234e-1/,p/1.7214206e-2/
      del=atan(tg*sin(p*(day-80.)))
      return
      end
