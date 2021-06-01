      subroutine conduc(l,nc,nl,na,alt,nh,dtets,par,stt,sft,sff,kpar)
      dimension alt(nh),par(kpar,na,nc),stt(nl,na,nc),sft(nl,na,nc),
     *sff(nl,na,nc)
      data e/1.60219e-20/,oe/1.76e7/,oi/3.11e2/,g10/-.30356/,
     *re/6371.02e5/,ci1/4.23e-10/,ci2/4.28e-10/,ci3/2.58e-10/,
     *ce1/1.82e-10/,ce11/3.6e-2/,ce2/2.33e-11/,ce21/1.21e-4/,
     *ce3/2.8e-10/,pi/3.14159265359/
      ie=(nc+1)/2
      do2i=1,nc
        tet=(i-1)*dtets/180.*pi
        ct=cos(tet)
        st=sin(tet)
        sk=sqrt(1.+3.*ct*ct)
        si=(ct+ct)/sk
        ci=st/sk
        sis=si*si
        cis=ci*ci
        do2j=1,na
          roq=(re/(re+alt(j)))**3
          b=-g10*roq*sk
          ee=e/b
          ome=oe*b
          omi=oi*b
          cm=par(4,j,i)
          co2=par(1,j,i)
          cn2=par(2,j,i)
          co=par(3,j,i)
c         tn=par(5,j,i)
          te=par(5,j,i)
c         b=sqrt(tn)
          b=sqrt(te)
          fi=ci1*co2+ci2*cn2+ci3*co
c         fe=ce1*(1.+ce11*b)*b*co2+ce2*(1.-ce21*tn)*tn*cn2+ce3*b*co
          fe=ce1*(1.+ce11*b)*b*co2+ce2*(1.-ce21*te)*te*cn2+ce3*b*co
          bu=ome*ome
          cu=omi*omi
          co2=fi*fi
          cn2=fe*fe
          b=1./(bu+cn2)
          co=1./(cu+co2)
          sp=cm*ee*(omi*fi*co+ome*fe*b)
          sh=cm*ee*(cu*co-bu*b)
          s0=cm*ee*(omi/fi+ome/fe)
          if(i.ne.ie)goto1
            stt(l,j,i)=s0
            sft(l,j,i)=0.
            sff(l,j,i)=sp+sh*sh/sp
            goto2
    1     continue
          b=1./(s0*sis+sp*cis)
          stt(l,j,i)=s0*sp*b
          sft(l,j,i)=s0*sh*si*b
          sff(l,j,i)=(s0*sp*sis+(sp*sp+sh*sh)*cis)*b
    2 continue
      return
      end

