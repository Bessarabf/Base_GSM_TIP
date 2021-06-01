      function  dela(beta,sk00r,sk00s)
       common/bldela/ap1,ap2,ap3,ap4,ap12av,ap36av,utsec
      g(ap)=(ap-4.)+(sk00r-1.)*
     *(ap-4.+(exp(-sk00s*(ap-4.))-1.)/sk00s)
      e=exp(-1.08e4*beta)
      e1=e**(amod(utsec,1.08e4)/1.08e4-0.5)
      s=1.+(((1.-e**19)*e1)/(1.-e))
      dela=g(ap1)+
     *   (g(ap2)*e+g(ap3)*e**2+g(ap4)*e**3+
     *   (g(ap12av)*e**4+g(ap36av)*e**12)*
     *   ((1.-e**8)/(1.-e))*e1)/s
      return
      end
