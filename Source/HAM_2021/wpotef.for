      subroutine wpotef(readfl,potef,n,kdf,kdu,
     *                   ldor,isp)
      dimension potef(n)
	allocatable pole(:)

      dimension kdf(20),kdu(20)
      logical readfl

	allocate (pole(ldor/4))

      mdor=ldor/4
      nfile=4
      nf=5
      if(nfile.le.8)nf=4
      isp=kdf(nfile)+1
      jj=kdu(4)
      jm=jj-1
      l=jj*mdor
      lm=l-mdor
      ll=n-lm
      if(readfl) go to 4
        k=1
        do 2 i=1,jm
          do 1 j=1,mdor
            pole(j)=potef(k)
            k=k+1
    1     continue
          write(nf,rec=isp)pole
          isp=isp+1
    2   continue
        do 3 j=1,ll
          pole(j)=potef(k)
          k=k+1
    3   continue
        write(nf,rec=isp)pole
        go to 8
    4 continue
        k=1
        do 6 i=1,jm
          read(nf,rec=isp)pole
          do 5 j=1,mdor
            potef(k)=pole(j)
            k=k+1
    5     continue
          isp=isp+1
    6   continue
        read(nf,rec=isp)pole
        do 7 j=1,ll
          potef(k)=pole(j)
          k=k+1
    7   continue
    8 continue
      deallocate (pole)
      return
      end
