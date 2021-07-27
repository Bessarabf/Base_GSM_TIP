      subroutine pecond(nl,nc,dtets,dfis,gamma,delta,psi)
      dimension gamma(nl,nc),delta(nl,nc),psi(nl,nc)
	allocatable a(:)
      allocate (a(nc))
	do 1 j=1,nc
        a(j)=90.-(j-1)*dtets
    1 continue
      d o 7i=1,nl
        fi=(i-1)*dfis
        print2,fi
    2   format(' ',53x,'long=',f5.0,' deg')
        j1=nc/10
        if(j1.ne.0)goto3
          i1=1
          i2=nc
          goto5
    3   continue
        do4j=1,j1
          i1=(j-1)*10
          i2=i1+10
          i1=i1+1
          print8,(a(k),k=i1,i2)
          print9,(gamma(i,k),k=i1,i2)
          print10,(delta(i,k),k=i1,i2)
          print11,(psi(i,k),k=i1,i2)
          print12
    4   continue
        j1=nc-j1*10
        if(j1.eq.0)goto6
          i1=nc-j1+1
          i2=nc
    5   continue
        print8,(a(k),k=i1,i2)
        print9,(gamma(i,k),k=i1,i2)
        print10,(delta(i,k),k=i1,i2)
        print11,(psi(i,k),k=i1,i2)
        print12
    6   continue
        print12
    7 continue
    8 format(' ','lat deg',3x,10f11.0)
    9 format(' ','stt mho',3x,1p10e11.2)
   10 format(' ','sff mho',3x,1p10e11.2)
   11 format(' ','sft mho',3x,1p10e11.2)
   12 format(' ',' ')

      deallocate(a)
      return
      end

