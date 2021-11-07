      subroutine libmanmf(nc,nl,na,alt,nh,dtets,dfis,alfa,
     *           beta,gamma,delta,psi,pef0,om,tet1,tet2,tet3,
     *           ae,be,ge,pef,u,nui,nvg)
      dimension alt(nh),alfa(nl,nc),beta(nl,nc),gamma(nl,nc),ae(nl),
     *          delta(nl,nc),psi(nl,nc),pef0(nl,nc),pef(nl,nc),be(nl),
     *          ge(nl),u(nui)
      data pi/3.14159265359/,re/6371.02e5/,pmi0/14.9/
      df=dfis/180.*pi
      dt=dtets/180.*pi
      dt1=dt*.5
      df=1./df
      dt=1./dt
      k=nc-1
      dfs=df*df*.5
      dts=dt*dt*.5
      df=df*.5
      dt=dt*.5
      m=(nc+1)/2
      mm=m-3
      mp=m+3
c      s1=(re+alt(na))/re*2.
      s1=(re+alt(na))/(re*2.)
      s2=s1*s1
      do4i=1,nl
        ip=i+1
        im=i-1
        if(i.eq.1)im=nl
        if(i.eq.nl)ip=1
        do3j=2,k
          if(j.eq.m)then
            call libme(i,ip,im,j,u,nui,nl,ae,be,ge,nc,pef0,pef,om,df)
c           pef(i,j)=0.
c           pef(i,j)=pef(i,j-1)
          else
            p=0.
            f=0.
            if(j.gt.mm.and.j.lt.mp)then
              l=nc-j+1
              n=j-1
              if(j.lt.m)then
                call libmu(i,ip,im,j,l,n,nl,nc,alfa,beta,gamma,delta,
     *          psi,u,nui,df,dfs,pef0,pef,om,p,f)
              else
                pef(i,j)=pef(i,l)
              end if
              goto3
            else
              if(j.le.mm)then
                tetn=(j-1)*dtets/180.*pi
              else
                tetn=(j-5)*dtets/180.*pi
              end if
              st=sin(tetn)
              sts=st*st
              pmi=(re+alt(na))/(re*sts)
              if(pmi.ge.pmi0)goto1
                l=nc-j+1
                if(j.gt.m)then
                  pef(i,j)=pef(i,l)
                  goto3
                else
c                  if(l.le.mm)tets=(l-1)*dtets/180.*pi
c                  if(l.ge.mp)tets=(l-5)*dtets/180.*pi
                  tets=(l-5)*dtets/180.*pi
                  call libmsn(i,ip,im,l,alfa,beta,gamma,delta,psi,dt,
     *            dts,df,dfs,pef0,nl,nc,p,f,s1,s2,dt1,tets)
c    *            dts,df,dfs,pef,nl,nc,p,f,s1,s2,dt1,tets)
                end if
    1         continue
              if(nvg.eq.0)goto2
                if(j.eq.4.or.j.eq.38)goto3
    2         continue
              call libmsn(i,ip,im,j,alfa,beta,gamma,delta,psi,dt,
     *        dts,df,dfs,pef0,nl,nc,p,f,s1,s2,dt1,tetn)
c    *        dts,df,dfs,pef,nl,nc,p,f,s1,s2,dt1,tetn)
              pef(i,j)=om/p*f+(1.-om)*pef0(i,j)
            end if
          end if
    3   continue
    4 continue
      s1=0.
      s2=0.
      do5i=1,nl
        s1=s1+pef(i,2)
        s2=s2+pef(i,k)
    5 continue
      s1=s1/nl
      s2=s2/nl
      do6i=1,nl
        pef(i,1)=s1
        pef(i,nc)=s2
    6 continue
      return
      end