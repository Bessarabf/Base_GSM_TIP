      subroutine psieq(l,nc,nl,na,alt,dtets,par,stt,sft,sff,kpar,psi,
     *              psie)
      dimension alt(na),par(kpar,na,nc),stt(nl,na,nc),psi(nl,nc),
     *         sft(nl,na,nc),sff(nl,na,nc),psie(nl)
      data g10/-.30356/,pi/3.14159265359/,re/6371.02e5/
      i=(nc+1)/2
c      roq=(re/(re+alt(1)))**3
c      b=-g10*roq
      b=bdip(alt(1),pi*.5)
      r=re+alt(16)
      vf=par(8,1,i)
      vr=par(6,1,i)
      tt=stt(l,1,i)
      ff=sff(l,1,i)
      h=alt(1)
      s1=0.
      s2=0.
c         substorm
c      nm=na-4
c         substorm
c         matias firster
      nm=na-5
c         matias firster
c	nm=na-4
      if(l.eq.1)then
        print2,' источник динамо-поля на экваторе',
     *         ' рассчитывается до ',nm,'-ой высотной точки'
    2   format(a34,a19,i3,a18)
      end if
      do1j=2,nm
c        roq=(re/(re+alt(j)))**3
c        bu=-g10*roq
        bu=bdip(alt(j),pi*.5)
        vfu=par(8,j,i)
	vru=par(6,j,i)
        ttu=stt(l,j,i)
	ffu=sff(l,j,i)
	hu=alt(j)
        dh=(hu-h)*.5
        s1=s1+(b*tt*vf+bu*ttu*vfu)*dh*2.
	s2=s2+(b*ff*vr+bu*ffu*vru)*dh
        b=bu
        vf=vfu
	vr=vru
        tt=ttu
	ff=ffu
        h=hu
    1 continue
      psi(l,i)=r*s1
c	psie(l)=0.
      psie(l)=r*s2
      return
      end
