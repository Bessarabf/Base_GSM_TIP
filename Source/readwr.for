c
c     READWR.f     11.03.93 Naumowa
c    cyclt:readwr
c    cyclt:cyclt1:readwr
c    cyclt:cyclt1:sklad:readwr
      subroutine readwr (readfl,nfile,kdf,kdu,nkd,nrazm,param,
     *                  ldor,isp,nf)
      dimension param(nrazm),kdu(nkd),kdf(nkd)
      logical readfl
      print *,' readwr: gins --> disk'
      print *,' readwr:  nf=',nf,'   nfile=',nfile
      mdor=ldor/4
      nk=mdor*kdu(nfile)
      if(nk.lt.nrazm)go to 1
        isp=kdf(nfile)+1
c       nf=4
c       if(nfile.gt.8)nf=5
        ind=1
        indk=mdor
        if(nrazm.lt.mdor)go to 8
          ncycl=nrazm/mdor
          nost=nrazm-ncycl*mdor
          if(nost.gt.0)ncycl=ncycl+1
          go to 9
    8   continue
          indk=nrazm
          ncycl=1
    9   continue
        if(readfl)go to 2
          do 5 i=1,ncycl
            write(nf,rec=isp)(param(j),j=ind,indk)
            isp=isp+1
            ind=indk+1
            indk=indk+mdor
            if(indk.gt.nrazm)indk=nrazm
    5     continue
          go to 3
    2   continue
          do 6 i=1,ncycl
            read(nf,rec=isp)(param(j),j=ind,indk)
            isp=isp+1
            ind=indk+1
            indk=indk+mdor
            if(indk.gt.nrazm)indk=nrazm
    6     continue
    3   continue
        go to 4
    1 continue
        print 900,nrazm,nfile,kdu(nfile),ldor
  900 format(' p/p readwr:  ****** ERROR !! *****'/
     *       '   nrazm=',i5,'      kdu(',i3,') =',i10)
        stop
    4 continue
      return
      end
