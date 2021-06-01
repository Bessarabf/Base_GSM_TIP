      function ds(altm,alt,altp,tetm,tetp)
      data re/6371.02e5/
      ds=sqrt((altp-altm)**2+(alt+re)**2*(tetp-tetm)**2)
      return
      end
