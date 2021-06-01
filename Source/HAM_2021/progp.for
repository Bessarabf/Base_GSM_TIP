      subroutine progp(ppo,pfu,ap,bp,ep,cp,dpp,nm,b00,
     *                g00,ed,fnp1)
      dimension ppo(nm),pfu(nm),ap(nm),cp(nm),ep(nm),
     *          dpp(nm),ed(nm),bp(nm)
      allocatable pb(:),pg(:)
	allocate (pb(nm),pg(nm))

      pb(1)=b00
      pg(1)=g00
      nmm=nm-1
      do 1 i=1,nmm
        if(abs(ep(i+1)).gt.0.3)go to 2
        ab=ed(i+1)+(1.-pb(i)*ep(i+1))*ap(i+1)
        pb(i+1)=(ep(i+1)*pb(i)-1.)/ab
        pg(i+1)=ep(i+1)*(pg(i)-pb(i)*bp(i))/ab
        go to 1
    2 continue
        ab=dpp(i+1)+(cp(i+1)-pb(i))*ap(i+1)
        pb(i+1)=(pb(i)-cp(i+1))/ab
        pg(i+1)=(pg(i)-pb(i)*bp(i))/ab
    1 continue
      ppo(nm)=(1.+pb(nm)*ap(nm))*
     *(fnp1+bp(nm))-pg(nm)*ap(nm)
      pfu(nm)=pg(nm)-pb(nm)*(bp(nm)+fnp1)
      do 3 j=2,nm
        i=nm-j+1
        ppo(i)=(1.+pb(i)*ap(i))*(ppo(i+1)+bp(i))
     *  -pg(i)*ap(i)
        pfu(i)=pg(i)-pb(i)*(ppo(i+1)+bp(i))
    3 continue
c     print*,'progp end'
      deallocate (pb,pg)

      return
      end

