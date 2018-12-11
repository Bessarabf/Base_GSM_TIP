      subroutine cint(kk,nvar,nqi,ntsl,kparn,ncn,qq,nxi,park,ns,
     *           pari,ni,nl,alfa,beta,sfti,psf,pst,rads,nh,par,kpars,ncs
     *          ,sip,sih,sib)
      dimension msum(45),ntsl(nqi),qq(nxi,nqi),park(ns),pari(ni),
     *          rads(nh),par(kpars,nh,ncs),
     *          sip(nl,ncn),sih(nl,ncn),sib(nl,ncn),
     *          alfa(nl,ncn),beta(nl,ncn),sfti(nl,ncn)
     *         ,psf(nl,ncn),pst(nl,ncn)
      data pi/3.14159265359/,re/6371.02e5/
      pmc=pi/180.
      msum(1)=0
      do 1 n=2,nqi
        nm=n-1
        msum(n)=msum(nm)+ntsl(nm)
    1 continue
      do 12 n=1,nqi
        nsum=msum(n)*kparn
        nt=ntsl(n)
        j1=ncn-n
        j2=n+1
        ls=0
        if(n.eq.1)goto3
          do2k=2,n
            ls=ls+ntsl(n-1)*2
    2     continue
   3   continue
        l=ls+1
          i1=1
          i3=(nt+1)/2
    5   continue
        i2=i1+1
        h=park(l)
        r=re+h
        ro=re/r
        ros=ro*ro
        roq=ros*ro
        t=park(l+1)*pmc
        ct=cos(t)
        st=sin(t)
        sks=1.+3.*ct*ct
        sk=sqrt(sks)
        sts=st*st
        b=bdip(h,t)
        q=qq(i1,n)
        ll=nsum+kparn*(i1-1)+1
        cm=pari(ll+3)
        co2=pari(ll)
        cn2=pari(ll+1)
        co=pari(ll+2)
        tn=pari(ll+4)
        vnu=pari(ll+6)
        vnv=pari(ll+7)
        call conducn(b,cm,co2,cn2,co,tn,sp,sh)
        s1=0.
        s2=0.
        s3=0.
        s4=0.
        s5=0.
        s9=0.
	  s10=0.
	  s11=0.
        do7i=i2,i3
          l=l+2
          h=park(l)
          rp=re+h
          rop=re/rp
          rosp=rop*rop
          roqp=rosp*rop
          t=park(l+1)*pmc
          ctp=cos(t)
          stp=sin(t)
          sksp=1.+3.*ctp*ctp
          skp=sqrt(sksp)
          stsp=stp*stp
          bp=bdip(h,t)
          qp=qq(i,n)
          ll=nsum+kparn*(i-1)+1
          cmp=pari(ll+3)
          co2p=pari(ll)
          cn2p=pari(ll+1)
          cop=pari(ll+2)
          tnp=pari(ll+4)
          vnup=pari(ll+6)
          vnvp=pari(ll+7)
          call conducn(bp,cmp,co2p,cn2p,cop,tnp,spp,shp)
          dq=(qp-q)*.5
          s1=s1+(sp*r*sts/ro+spp*rp*stsp/rop)*dq
          s2=s2+(sp*r/(roq*sts*sks)+spp*rp/(roqp*stsp*sksp))*dq
          s3=s3+(sh*r/(ros*sk)+shp*rp/(rosp*skp))*dq
	    s9=s9+(sp/(roq*sk)+spp/(roqp*skp))*dq*re
	    s10=s10+(sh/(roq*sk)+shp/(roqp*skp))*dq*re
          if(nvar.eq.2)goto6
            s4=s4+((sp*vnu-sh*vnv)*b*r*r/(roq*st*sks)+(spp*vnup-
     *      shp*vnvp)*bp*rp*rp/(roqp*stp*sksp))*dq
            s5=s5+((sh*vnu+sp*vnv)*b*r*r*st/(ros*sk)+
     *      (shp*vnup+spp*vnvp)*bp*rp*rp*stp/(rosp*skp))*dq
	      s11=s11+((sp*vnu-sh*vnv)*b/(roq*sk)+
     *	           (spp*vnup-shp*vnvp)*bp/(roqp*skp))*dq*re
    6     continue
          ct=ctp
          st=stp
          sks=sksp
          sk=skp
          sts=stsp
          r=rp
          ro=rop
          ros=rosp
          roq=roqp
          b=bp
          q=qp
          vnu=vnup
          vnv=vnvp
	    sp=spp
	    sh=shp
    7   continue
        if(i1.ne.1)goto8
          alfa(kk,j1)=s1
          beta(kk,j1)=s2
          sfti(kk,j1)=s3
          psf(kk,j1)=s4
          pst(kk,j1)=s5
	    sip(kk,j1)=s9
	    sih(kk,j1)=s10
	    sib(kk,j1)=s11
          goto9
    8   continue
          alfa(kk,j2)=s1
          beta(kk,j2)=s2
          sfti(kk,j2)=s3
          psf(kk,j2)=s4
          pst(kk,j2)=s5
	    sip(kk,j2)=s9
	    sih(kk,j2)=s10
	    sib(kk,j2)=s11
    9   continue
        if(i3.eq.nt)goto12
          i1=(nt+2)/2
          if(i1.eq.i3)goto10
            l=l+2
   10     continue
          i3=nt
          goto5
   12 continue
      i=(ncn+1)/2
      pst(kk,i)=0.
      pst(kk,1)=0.
      pst(kk,ncn)=0.
      alfa(kk,i)=0.
      alfa(kk,1)=0.
      alfa(kk,ncn)=0.
      sfti(kk,i)=0.
	sip(kk,1)=0.
	sih(kk,1)=0.
	sib(kk,1)=0.
	sip(kk,ncn)=0.
	sih(kk,ncn)=0.
	sib(kk,ncn)=0.
      return
      end
