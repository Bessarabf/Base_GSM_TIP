      subroutine libmu(i,ip,im,j,l,n,nl,nc,alfa,beta,gamma,delta,psi,
     *u,nui,df,dfs,pef0,pef,om,p,f)
      dimension alfa(nl,nc),beta(nl,nc),gamma(nl,nc),u(nui),
     *delta(nl,nc),psi(nl,nc),pef0(nl,nc),pef(nl,nc)
      jp=j+1
      jm=j-1
      lp=l-1
      lm=l+1
      du=u(n+1)-u(n)
      du=1./du
      dus=du*du*.5
      du=du*.5
      e=alfa(i,l)
      g=(e+alfa(i,lp))*dus
      h=gamma(i,l)*du
      a=abs(h)
      b=g+h+a
      p=p+b
c     f=f+b*pef(i,lp)
      f=f+b*pef0(i,lp)
      g=(e+alfa(i,lm))*dus
      b=g-h+a
      p=p+b
c     f=f+b*pef(i,lm)
      f=f+b*pef0(i,lm)
      e=beta(i,l)
      h=delta(i,l)*df
      g=(e+beta(ip,l))*dfs
      a=abs(h)
      b=g-h+a
      p=p+b
c     f=f+b*pef(ip,l)
      f=f+b*pef0(ip,l)
      g=(e+beta(im,l))*dfs
      b=g+h+a
      p=p+b
c     f=f+b*pef(im,l)-psi(i,l)
      f=f+b*pef0(im,l)-psi(i,l)
      e=alfa(i,j)
      h=gamma(i,j)*du
      g=(e+alfa(i,jp))*dus
      a=abs(h)
      b=g+h+a
      p=p+b
c     f=f+b*pef(i,jp)
      f=f+b*pef0(i,jp)
      g=(e+alfa(i,jm))*dus
      b=g-h+a
      p=p+b
c     f=f+b*pef(i,jm)
      f=f+b*pef0(i,jm)
      e=beta(i,j)
      h=delta(i,j)*df
      g=(e+beta(ip,j))*dfs
      a=abs(h)
      b=g-h+a
      p=p+b
c     f=f+b*pef(ip,j)
      f=f+b*pef0(ip,j)
      g=(e+beta(im,j))*dfs
      b=g+h+a
      p=p+b
c     f=f+b*pef(im,j)
      f=f+b*pef0(im,j)
      f=f-psi(i,j)
      pef(i,j)=om/p*f+(1.-om)*pef0(i,j)
      return
      end
