      subroutine sftpol(ll,rads,nh,na,par,kpars,ncs,ncn,sfti,nl)
      dimension rads(nh),par(kpars,nh,ncs),sfti(nl,ncn)
      data re/6371.02e5/
      i=1
    1 continue
      h=rads(1)
      r=re+h
      ro=(re/r)**2
      t=0.
      b=bdip(h,t)
      cm=par(6,1,i)
      co2=par(1,1,i)
      cn2=par(2,1,i)
      co=par(3,1,i)
      tn=par(7,1,i)
      call conducn(b,cm,co2,cn2,co,tn,sp,sh)
      s=0.
      do2j=2,na
        h=rads(j)
        rp=re+h
        rop=(re/rp)**2
        bp=bdip(h,t)
        cmp=par(6,j,i)
        co2p=par(1,j,i)
        cn2p=par(2,j,i)
        cop=par(3,j,i)
        tnp=par(7,j,i)
        call conducn(bp,cmp,co2p,cn2p,cop,tnp,spp,shp)
        dq=abs(rop-ro)*.25
        s=s+(sh*r/ro+shp*rp/rop)*dq
        r=rp
        ro=rop
        b=bp
        cm=cmp
        co2=co2p
        cn2=cn2p
        co=cop
        tn=tnp
    2 continue
      if(i.eq.ncs)goto3
        sfti(ll,1)=s
        i=ncs
        goto1
    3 continue
      sfti(ll,ncn)=s
      return
      end

