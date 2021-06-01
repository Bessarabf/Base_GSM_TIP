      subroutine rpco(dut,hut,iput,keut,koob,kut,nara,naou,sut,tut,
     *                isat)
c   out: dut  - longitude of dot orbit
c        hut  - height of dot orbit
c        iput - number of dot height profile
c        keut - length of record
c        koob - number of object
c        kut  - number of dot orbit
c        nara - name of input file object
c        naou - name of output file object
c        sut  - latitude of dot orbit
c        tut  - UT of dot orbit
  
      dimension dut(isat,100),hut(isat,100),iput(100),keut(100),
     *          kut(100),sut(isat,100),tut(isat,100)
      character naob(100)*4,nara(100)*25,naou(100)*25
      open(9,file='fina',status='old')
      i=0
    1 continue
        i=i+1
        read(9,'(a)',end=2) naob(i)
        goto 1
    2 continue
!      koob=i-1    ! „«ï NDP
      koob=i-2    ! „«ï FL32
c     print *, kut
      do i=1,koob
        nara(i)='raf/raf'//naob(i)
        naou(i)='ouf'//naob(i)
c       print*,naob(i),nara(i),naou(i),koob
      end do
      close(9)
      do i=1,koob
        open(10,file=nara(i),status='old')
        read(10,'(i4)') kout
c       print*,nara(i),i,kout
        kut(i)=kout
        ha=0.
        do k=1,kout
          read(10,'(i4,f9.3,2f8.2,f10.2)') ko,tut(k,i),sut(k,i),
     *                                     dut(k,i), h
          if(h.gt.ha) ha=h
          hut(k,i)=h
        end do
        close(10)
        call camu(ha,ip,kec)
        iput(i)=ip
        keut(i)=kec
c       print*,ha,ip,kec
      end do
      do i=1,koob
        kout=kut(i)
        a=tut(1,i)
        do k=2,kout
          b=tut(k,i)
          if(b.lt.a)then
            tut(k,i)=b+24.
            a=tut(k,i)
          else
            a=b
          end if
        end do
      end do
      return
      end
*****************************************************
      subroutine camu(ha,ip,kec)
      dimension he(70),hei(70)
      data hei/ 175.3, 187.8, 201.6, 216.8, 233.5, 251.8, 272.0, 294.2,
     *          318.6, 345.5, 375.0, 407.5, 443.3, 482.6, 525.9, 573.5,
     *          625.8, 683.4, 746.8, 816.4, 893.1, 977.4,1070.1,1172.1,
     *         1284.3,1407.8,1543.6,1692.9,1857.2,2037.9,2236.7,2455.4,
     *         2695.9,2960.5,3251.6,3571.7,3923.9,4311.3,4737.4,5206.2,
     *         5721.8,6289.0,6912.9,7599.1,8354.1,9184.5,10097.9,
     *        11102.7,12208.0,13423.8,14761.1,16232.3,17850.5,19630.5,
     *        21588.6,23742.4,26111.7,28717.9,31584.6,34738.1,38206.9,
     *        42022.6,46219.9,50836.9,55915.5,61502.1,67647.3,74407.0,
     *        81842.8,88557.2/
      n=70
      do i=1,n
        he(i)=hei(i)
      end do
      ip=j52(n,ha,he)+2
      if(ip.gt.70)ip=70
      kec=(ip*8+480+9)*4
      return
      end
******************************************************
      function j52(n,u,x)
      dimension x(n)
      i=n
      s=x(i)
      i=n-1
      if(u.ge.s) goto 2
      i=1
      s=x(i+1)
      if(u.lt.s) goto 2
      j=n+1
    1 continue
        k=(i+j)/2
        s=x(k)
        if(u.lt.s) j=k
        if(u.ge.s) i=k
      if(j.gt.i+1) goto 1
    2 continue
      j52=i
      return
      end
