      subroutine libmsn(i,ip,im,j,alfa,beta,gamma,delta,psi,
     *dt,dts,df,dfs,pef0,nl,nc,p,f,s1,s2,dt1,tet)
      dimension alfa(nl,nc),beta(nl,nc),gamma(nl,nc),
     *delta(nl,nc),psi(nl,nc),pef0(nl,nc)
      jp=j+1
      jm=j-1
      e=alfa(i,j)
      g=(e+alfa(i,jp))*dts
      s3=sin(tet)*cos(tet)
c      ct=cos(tet)
c      s3=(1.+3.*ct*ct)*sin(tet)/ct
      tets=tet+dt1
      s4=sin(tets)*cos(tets)
c      ct=cos(tets)
c      s4=(1.+3.*ct*ct)*sin(tets)/ct
      s5=s2/s3
      g=g*s5/s4
      h=gamma(i,j)*dt*s1/s3
c      h=gamma(i,j)*dt*s1/s4
      a=abs(h)
      b=g+h+a
      p=p+b
      f=f+b*pef0(i,jp)
      g=(e+alfa(i,jm))*dts
      tets=tet-dt1
      s4=sin(tets)*cos(tets)
c      ct=cos(tets)
c      s4=(1.+3.*ct*ct)*sin(tets)/ct
      g=g*s5/s4
      b=g-h+a
      p=p+b
      f=f+b*pef0(i,jm)
      e=beta(i,j)
      h=delta(i,j)*df
      g=(e+beta(ip,j))*dfs
      a=abs(h)
      b=g-h+a
      p=p+b
      f=f+b*pef0(ip,j)
      g=(e+beta(im,j))*dfs
      b=g+h+a
      p=p+b
      f=f+b*pef0(im,j)
      f=f-psi(i,j)
      return
      end
