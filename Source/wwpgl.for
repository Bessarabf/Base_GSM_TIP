c                WWPGL  (cyclt1:potnew:wwpgl )
c             -----------------------------------
      subroutine wwpgl (pglo,kpars,nh,its,ids,par,dolg,ddolgs)
c                vhod - pglo,kpars,nh,its,ids,par,dolg,ddolgs
c                vihod - par
      dimension pglo(kpars,nh,its,ids),par(kpars,nh,its)
      nn=dolg
      nd=ddolgs
      if(nn.ge.360)nn=nn-360
      i=nn/nd+1
      do 2 j=1,its
        do 3 k=1,nh
          do 4 n=1,kpars
             par(n,k,j)=pglo(n,k,j,i)
    4     continue
    3   continue
    2 continue
      return
      end
