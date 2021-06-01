      subroutine degaps(nvar,nc,nl,na,alt,nh,dtets,dfis,gamma,
     *delta,psi,pst,psf,sfti)
      dimension alt(nh),gamma(nl,nc),delta(nl,nc),psi(nl,nc),
     *pst(nl,nc),psf(nl,nc),sfti(nl,nc)
      data pi/3.14159265359/,re/6371.02e5/
      df=(dfis+dfis)/180.*pi
      dt=(dtets+dtets)/180.*pi
      df=1./df
      dt=1./dt
      r=re+alt(na)
      k=nc-1
      k1=(nc+1)/2
      do2i=2,k
        ip=i+1
        im=i-1
        do2j=1,nl
          jp=j+1
          jm=j-1
          if(j.eq.1)jm=nl
          if(j.eq.nl)jp=1
          gamma(j,i)=(sfti(jp,i)-sfti(jm,i))*df
          delta(j,i)=(sfti(j,ip)-sfti(j,im))*dt
          if(nvar.eq.2)goto1
            if(i.eq.k1)goto2
              psi(j,i)=r*((pst(j,ip)-pst(j,im))*dt+(psf(jp,i)-
     *        psf(jm,i))*df)
              goto2
    1     continue
          psi(j,i)=0.
    2 continue
      return
      end

