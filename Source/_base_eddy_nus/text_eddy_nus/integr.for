      subroutine integr(l,nc,nl,na,alt,nh,dtets,par,stt,sft,sff,
     *           alfa,beta,sfti,pst,psf,kpar,gamma,delta,psi,
     *           sid,sift,siff,sitt)
      dimension alt(nh),par(kpar,na,nc),stt(nl,na,nc),sft(nl,na,nc),
     *          sff(nl,na,nc),alfa(nl,nc),beta(nl,nc),sfti(nl,nc),
     *          pst(nl,nc),psf(nl,nc),gamma(nl,nc),delta(nl,nc),
     *          psi(nl,nc),sid(nl,nc),siff(nl,nc),sift(nl,nc),
     *          sitt(nl,nc)
      data pi/3.14159265359/,re/6371.02e5/
      pmc=pi/180.
      do4i=1,nc
        tet=(i-1)*dtets*pmc
        ct=cos(tet)
        st=sin(tet)
        sk=sqrt(1.+3.*ct*ct)
        si=(ct+ct)/sk
        ci=st/sk
	  b=bdip(alt(1),tet)
        vt=par(7,1,i)
        vr=par(6,1,i)
        vf=par(8,1,i)
        tt=stt(l,1,i)
        ft=sft(l,1,i)
        h=alt(1)
        s1=0.
        s3=0.
        s5=0.
        s2=0.
        if(i.eq.1.or.i.eq.nc)goto1
          s4=0.
    1   continue
        ff=sff(l,1,i)
        do3j=2,na
          bu=bdip(alt(j),tet)
          vtu=par(7,j,i)
          vru=par(6,j,i)
          vfu=par(8,j,i)
          ttu=stt(l,j,i)
          ftu=sft(l,j,i)
          hu=alt(j)
          dh=(hu-h)*.5
          s1=s1+(tt+ttu)*dh
          c1=-(b*tt*vf+bu*ttu*vfu)*dh*si
          c2=-(b*ft*vt+bu*ftu*vtu)*dh*si
          c3=(b*ft*vr+bu*ftu*vru)*dh*ci
          s3=s3+c1+c2+c3
          s5=s5+(ft+ftu)*dh
          ffu=sff(l,j,i)
          s2=s2+(ff+ffu)*dh
          if(i.eq.1.or.i.eq.nc)goto2
            c4=-(b*ft*vf+bu*ftu*vfu)*dh*si
            c5=(b*ff*vt+bu*ffu*vtu)*dh*si
            c6=-(b*ff*vr+bu*ffu*vru)*dh*ci
            s4=s4+c4+c5+c6
    2     continue
          ff=ffu
          b=bu
          vt=vtu
          vr=vru
          vf=vfu
          tt=ttu
          ft=ftu
          h=hu
    3   continue
        gamma(l,i)=s1*1.e9
	  sitt(l,i)=s1
        alfa(l,i)=s1*st
        pst(l,i)=s3*st
        sfti(l,i)=s5
        delta(l,i)=s2*1.e9
	  sift(l,i)=s5
	  siff(l,i)=s2
        if(i.eq.1.or.i.eq.nc)goto4
          beta(l,i)=s2/st
          psf(l,i)=s4
	    sid(l,i)=s4
    4 continue
      return
      end
