      subroutine zamion(gins,ins,its,ids,nh,uts,mass)

      dimension gins(ins,nh,its,ids),mass(30)
	allocatable s(:)
	allocate (s(ids))
      if(mass(13).ne.0)go to 1
c     . . .
        iutss=nint(uts/3600.)
c       iutss=uts/3600.
c     . . .
   11   if(iutss.lt.24) go to 10
          iutss=iutss-24
          go to 11
   10   continue
        nd=iutss-mass(14)
        do 2 np=1,ins
          do 3 k=1,nh
            do 4 i=1,its
              do 5 j=1,ids
                jn=j+nd
                if(jn.gt.ids)jn=jn-ids
                if(jn.lt.1)jn=jn+ids
                s(j)=gins(np,k,i,jn)
    5         continue
              do 6 j=1,ids
                gins(np,k,i,j)=s(j)
    6         continue
    4       continue
    3     continue
    2   continue
    1 continue
      deallocate (s)

      return
      end
                                                                                        
