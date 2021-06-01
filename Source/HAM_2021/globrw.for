c  version GSM  15.04.14
       subroutine globrw(nfile,readfl,pgl1,pole,kpar,ldor,
     *                  nh,kdf,isp,npgl,its,ids,mass)
 
      dimension pgl1(npgl),kdf(20),mass(30),pole(ldor/4)
      logical readfl
c
  667 format(' globrw - begin')
  767 format(' globrw - end')
  700 format(' globrw :   ********  error  ********'/
     *   '  npg=',i8,'  >   npgl=',i8,'  !!!!!!  STOP  !')
c
      print 667

      npg=kpar*nh*its*ids
      if(npgl.lt.npg) go to 9
        nf=5
        if(nfile.le.8) nf=4
        isp=kdf(nfile)+1
        mdor=ldor/4
        ndor=npgl/mdor
        nost=npgl-ndor*mdor
        if(readfl)go to 4
          k=1
          do 2 j=1,ndor
            do 1 i=1,mdor
              pole(i)=pgl1(k)
              k=k+1
    1       continue
            write(nf,rec=isp)pole
            isp=isp+1
    2     continue
          read(nf,rec=isp)pole
          do 3 i=1,nost
            pole(i)=pgl1(k)
            k=k+1
    3     continue
          write(nf,rec=isp)pole
          go to 8
    4   continue
          k=1
          do 6 j=1,ndor
            read(nf,rec=isp) pole
            do 5 i=1,mdor
              pgl1(k)=pole(i)
              k=k+1
    5       continue
            isp=isp+1
    6     continue
          if(nost.eq.0)go to 11
            read(nf,rec=isp)pole
            do 7 i=1,nost
              pgl1(k)=pole(i)
              k=k+1
    7       continue
   11     continue
    8   continue
        go to 10
    9 continue
        print 700,npg,npgl
        stop
   10 continue
      print 767

      return
      end
 
