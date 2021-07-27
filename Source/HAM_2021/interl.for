      subroutine interl(yz,z,nz,yx,x,nx)
c
c     . . . subroutine calculated array  yz on
c     argument array z
c     elments array z stay on increase order
c
      dimension yz(nz),z(nz),x(nx),yx(nx)
      data key/0/
      key=key+1
c     . . . region definition of arguments
      kx=1
      kz=1
      if(z(1).lt.x(1).or.z(nz).gt.x(nx)) then
        print 100,z(1),x(1),z(nz),x(nx)
        return
      endif
    1 kx=kx+1
    2 if(.not.(z(kz).ge.x(kx-1).and.
     *         z(kz).le.x(kx))) go to 3
       yz(kz)=(yx(kx)-yx(kx-1))/(x(kx)-x(kx-1))*
     *        (z(kz)-x(kx-1))+yx(kx-1)
       kz=kz+1
       if(kz.gt.nz)
     * return
       go to 2
    3 continue
      if(kx.eq.nx)
     * return
      go to 1
  100 format(' interval not found',4(1pe10.3))
      return
      end

