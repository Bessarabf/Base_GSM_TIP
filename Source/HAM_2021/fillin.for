      subroutine fillin(ne,j,m,ntsl,nl,cio,cih,cihe,vio,vih,vihe,
     *           ti,te,vdv,vdu,par,nr)
      dimension ntsl(nl),cio(*),cih(*),cihe(*),vio(*),
     *          vih(*),vihe(*),ti(*),te(*),vdv(*),vdu(*),
     *          par(nr)
      i1=0
      if(j.eq.1)goto2
        i2=j-1
        do1i=1,i2
          i1=i1+ntsl(i)
    1   continue
    2 continue
      n=ntsl(j)
      do 5 i=1,n
        k=(i1+i-1)*m+1
        goto(3,4),ne
    3   continue
          par(k)=cio(i)
          par(k+1)=cih(i)
          par(k+2)=cihe(i)
          par(k+3)=vio(i)
          par(k+4)=vih(i)
          par(k+5)=vihe(i)
          par(k+6)=ti(i)
          par(k+7)=te(i)
          goto5
    4   continue
          par(k)=vdv(i)
          par(k+1)=vdu(i)
    5 continue
      return
      end

