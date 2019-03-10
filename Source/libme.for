      subroutine libme(i,ip,im,j,u,nui,nl,ae,be,ge,nc,pef0,pef,om,df)
      dimension u(nui),ae(nl),be(nl),ge(nl),pef0(nl,nc),pef(nl,nc)
      jm=j-1
      a=1./(u(nui)-u(nui-1))
      b=-be(i)*a
      e=ae(i)
      a=abs(e)
      g=a+e
      e=e-a
c     pef(i,j)=1./(om*(b+(a+a)*df)+1.)*((b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))*om+pef0(i,j))
c     pef(i,j)=om/(om*(b+(a+a)*df)+4.)*((b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))*2.+(4.-om*(b+(a+a)*df))/om*pef0(i,j))
c     pef(i,j)=(.5*om/(b+(a+a)*df)*(b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))+(1.-om*.25)*pef0(i,j))/(1.+.25*om)
c     pef(i,j)=om/(b+(a+a)*df)*(b*pef(i,jm)+(g*pef(ip,j)-
c    *e*pef(im,j))*df-ge(i))+(1.-om)*pef0(i,j)
      pef(i,j)=om/(b+(a+a)*df)*(b*pef0(i,jm)+(g*pef0(ip,j)-
     *e*pef0(im,j))*df-ge(i))+(1.-om)*pef0(i,j)
      return
      end
