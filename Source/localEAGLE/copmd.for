      subroutine copmd(nfr,nfw,kdf,isp,ldor,kdu,par,nr,mass,mast)
      dimension par(nr),mast(40),mass(30),kdu(20),kdf(20)
  900 format(' copmd:',20i5)
  901 format(' copmd:  ********  error *********'/
     *       ' размеры участков на дисках',
     *       ' не совпадают'/
     *  ' kdu(',i2,')=',i7,'   kdu(',i2,')=',i7,' !!!!!!! STOP !')
  902 format(' copmd:  ********  error ********** '/
     *       ' размерность области I/O меньше mdor'/
     *         ' nr = ',i10,'  <   mdor=',i5,'  !!!!!! STOP  !')
  903 format(' p/p copmd: nfr=',i5,'  ---> ','  nfw=',i5)
      
	
	mdor=ldor/4
      if(nfr.eq.14.and.mast(13).eq.0)go to 9
        if(nfw.eq.14.and.mast(13).eq.0)nfw=13
        
        print 900,nfr,nfw
        if(nr.lt.mdor) go to 3
          i1=kdf(nfr)+1
      	i2=kdf(nfw)+1
          kd1=kdu(nfr)
          kd2=kdu(nfw)
          if(kd1.ne.kd2) go to 2
            nf1=5
            nf2=5
            if(nfr.le.8) nf1=4
            if(nfw.le.8) nf2=4
            do 1 j=1,kd1
              read(nf1,rec=i1)(par (i),i=1,mdor)
              write(nf2,rec=i2)(par (i),i=1,mdor)
              i1=i1+1
              i2=i2+1
    1       continue
            print 903,nfr,nfw
            go to 5
    2     continue
            print 901,nfr,kd1,nfw,kd2
            stop
    5     continue
          go to 8
    3   continue
          print 902,nr,mdor
          stop
    8   continue
        print *,' copmd - end'
    9 continue
      return
      end

