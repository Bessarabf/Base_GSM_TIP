c
      subroutine gstrf0(pgl,rads,kpars,nh,its,ids
     *                 ,dtets,ddolg,m usl)
      dimension pgl(kpars,nh,its,ids)
     *         ,rads(nh)
      data om,pi/7.27e-5,3.1416/,bk/1.38e-16/,
     *     am1,am2,am3/53.12e-24,46.56e-24,26.56e-24/
cc      Ветер на нижней границе =0 всегда!!!
c     if(.not.(m usl.eq.0)) go to 21
       do 22 i=1,its
        do 23 j=1,ids
         pgl(11,1,i,j)=0.
         pgl(12,1,i,j)=0.
   23   continue
   22  continue
       return
c  21 continue
c     ns=its-1
c     nsp=ns/2+1
c     nsp1=nsp-3
c     nsp2=nsp+3
c     do 111 j = 1 , ids
c      j1=j-1
c      j2=j+1
c      dolg=(j-1)*ddolg
c      if(j.eq.1) j1=ids
c      if(j.eq.ids) j2=1
c      dolgrd=dolg*pi/180.
c      dtetsr=dtets*pi/180.
c      ddolgr=ddolg*pi/180.
c      r=rads(1)+6.371e+8
c      do 3 i = 1 , its
c       pgl(11,1,i,j)=0.
c       pgl(12,1,i,j)=0.
c 3    continue
c      do 1 i = 2 , ns
c       tet=(i-1)*dtets*pi/180.
c       ro=am1*pgl(1,1,i,j)+am2*pgl(2,1,i,j)+am3*pgl(3,1,i,j)
c     di=(pgl(1,1,i+1,j)+pgl(2,1,i+1,j)+pgl(3,1,i+1,j))*pgl(7,1,i+1,j)
c     di1=(pgl(1,1,i-1,j)+pgl(2,1,i-1,j)+pgl(3,1,i-1,j))*pgl(7,1,i-1,j)
c     dj=(pgl(1,1,i,j1)+pgl(2,1,i,j1)+pgl(3,1,i,j1))*pgl(7,1,i,j1)
c     dj1=(pgl(1,1,i,j2)+pgl(2,1,i,j2)+pgl(3,1,i,j2))*pgl(7,1,i,j2)
c      csg=0.9807*cos(tet)-0.1956*sin(tet)*cos(dolgrd)
c      x=2.*om*csg
c      if(i.lt.nsp1.or.i.gt.nsp2) go to 11
c      if(abs(csg).lt.0.2) go to 1
c 11   pgl(11,1,i,j)=-bk*(dj1-dj)/x/r/sin(tet)/ddolgr/2./ro
c      pgl(12,1,i,j)=bk*(di-di1)/x/r/dtetsr/2./ro
c 1   continue
c     pgl(12,1,1,j)=2.*pgl(12,1,2,j)-pgl(12,1,3,j)
c     pgl(12,1,its,j)=2.*pgl(12,1,ns,j)-pgl(12,1,ns-1,j)
c     pgl(11,1,1,j)=2.*pgl(11,1,2,j)-pgl(11,1,3,j)
c     pgl(11,1,its,j)=2.*pgl(11,1,ns,j)-pgl(11,1,ns-1,j)
c        do 5 i = 1 , its
c       if(pgl(11,1,i,j).ne.0.) go to 5
c       pgl(11,1,i,j)=(pgl(11,1,i+1,j)+pgl(11,1,i-1,j))/2.
c       if(pgl(12,1,i,j).ne.0.) go to 5
c       pgl(12,1,i,j)=(pgl(12,1,i+1,j)+pgl(12,1,i-1,j))/2.
c 5    continue
c111   continue
c      return
       end
