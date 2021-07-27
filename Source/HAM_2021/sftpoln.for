      subroutine sftpoln(ll,rads,nh,na,par,kpars,ncs,ncn,sfti,nl)
      dimension rads(nh),par(kpars,nh,ncs),sfti(nl,ncn)
      sfti(ll,1)=0.
      sfti(ll,ncn)=0.
      return
      end
