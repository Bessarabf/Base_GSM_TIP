      subroutine intpar(ncs,nc,dtets,dtet,par,kpars,nh,par1,kpar,na)
      dimension par(kpars,nh,ncs),par1(kpar,na,nc)
      j1=1
      j2=2
      tet1=0.
      tet2=dtets
      j3=nc-1
      do4j=1,j3
        te=(j-1)*dtet
    1   continue
        if(te.ge.tet1.and.te.lt.tet2)goto2
          tet1=tet2
          tet2=tet2+dtets
          j1=j2
          j2=j2+1
          goto1
    2   continue
        dt=(te-tet1)/dtets
        do3ip=1,kpar
          if(ip.lt.4)ips=ip
c                PC
          if(ip.eq.4.or.ip.eq.5)ips=ip+2
          if(ip.gt.5)ips=ip+4
c                PC
c                Labtam
c         if(ip.eq.4)ips=ip+2
c         if(ip.ge.5)ips=ip+4
c                Labtam
          do3k=1,na
            df=par(ips,k,j2)-par(ips,k,j1)
            par1(ip,k,j)=par(ips,k,j1)+df*dt
    3   continue
    4 continue
      do5ip=1,kpar
        if(ip.lt.4)ips=ip
c                PC
        if(ip.eq.4.or.ip.eq.5)ips=ip+2
        if(ip.gt.5)ips=ip+4
c                PC
c                Labtam
c       if(ip.eq.4)ips=ip+2
c       if(ip.ge.5)ips=ip+4
c                Labtam
        do5k=1,na
          par1(ip,k,nc)=par(ips,k,ncs)
    5 continue
      return
      end

