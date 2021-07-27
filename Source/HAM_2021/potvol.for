c . . . ver. 2012 
      subroutine potvol(delta,ut,kp,bmpy,alt,ddolg,dtet,potef,
     *                 idt,nl2,ntr,b,c,isp,readfl,kdf,kdu,ldor)
      double precision pi
      real kp
      dimension potef(ntr,idt,nl2),kdf(20),kdu(20)
      logical readfl
      data pi/3.14159265359d0/
      i2=idt/2
      print 7
      j1=(nl2+1)/2
      do 6 i=1,idt
        fi=(i-1)*pi/i2
        do 6 j=1,nl2
          if(j.ne.1)go to 1
            tet=0.
            go to 5
    1     continue
          if(j.ne.nl2)go to 2
            tet=pi
            go to 5
    2     continue
          if(j.ne.j1)go to 3
            tet=pi*0.5
            go to 5
    3     continue
          if(j.gt.j1)go to 4
            tet=(c-(j-j1+1)*dtet)*pi/180.
            go to 5
    4     continue
            tet=(b-(j-2)*dtet)*pi/180.
    5  continue
        potef(16,i,j)=pefvol(delta,ut,kp,bmpy,alt,tet,fi,j,j1,nl2)
    6 continue
    7 format(' potvol')
      na=16
      readfl=.false.
      nanl=na*idt*nl2
      call wpotef(readfl,potef,nanl,kdf,kdu,ldor,isp)
      
      pmaxn=potef(na,1,1)
      pminn=potef(na,1,1)
      pmaxs=potef(na,1,nl2)
      pmins=potef(na,1,nl2)
      do8j=1,j1
        k=nl2-j+1
        do8i=1,idt
          pn=potef(na,i,j)
          ps=potef(na,i,k)
          if(pmaxn.lt.pn)pmaxn=pn
          if(pminn.gt.pn)pminn=pn
          if(pmaxs.lt.ps)pmaxs=ps
          if(pmins.gt.ps)pmins=ps
    8 continue
      dpn=(pmaxn-pminn)*1.e-11
      dps=(pmaxs-pmins)*1.e-11
      print  9,dps,dpn
    9 format(' ',35x,'dps=',1pe9.2,' kV',10x,'dpn=',1pe9.2,' kV')
      return
      end
