      subroutine filpef(idol,j,ntsl,ncs,pef,nl,nc,pefgl,nr)
      dimension ntsl(ncs),pefgl(nr),pef(nl,nc)
      i1=0
      if(j.eq.1)goto1
        i2=j-1
        do i=1,i2
          i1=i1+ntsl(i)
        end do
    1 continue
      n=ntsl(j)
      do i=1,n
        k=i1+i
        if(j.lt.5)then
          n1=n/2
          if(i.le.n1)then
            j1=nc-j+1
          else
            j1=j
          end if
        else
          j1=j
        end if
        pefgl(k)=pef(idol,j1)
      end do
      return
      end
