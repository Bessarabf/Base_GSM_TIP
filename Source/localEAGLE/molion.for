





      subroutine molion(par,kpars,nh,its,pari,ins,mass,dts,ntr)
      dimension par(kpars,nh,its),mass(30),pari(ins,nh,its)
      real    l1,l2,moli,lamo
      data alfa/4.2e-7/,l1/1.6e-11/,l2/6.e-13/
      do 1 ig=1,its
        do 2 i=1,nh
          moli=par(6,i,ig)
      if(par(9,i,ig).le.0.1)print 900,par(9,i,ig),i,ig
  900 format(' par(9,i,ig)=',g12.4,2i5)
          alt=alfa*(300./par(9,i,ig))
          qs=par(13,i,ig)+par(14,i,ig)+par(15,i,ig)
          if (i.ge.ntr) go to 5
            qs=qs+par(16,i,ig)
            go to 6
    5     continue
          vrn=par(10,i,ig)
          vtn=par(11,i,ig)
          vln=par(12,i,ig)
          vior=pari(2,i,ig)
          viot=pari(3,i,ig)
          viol=pari(4,i,ig)
          tn=par(7,i,ig)
c
          tim=tn
          if(i.ge.16)tim=pari(6,i,ig)
          l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
          l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,tn)
c         te=par(9,i,ig)
c         l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c         l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol,te)
c         l1=lamo(1,tn,vrn,vtn,vln,tim,vior,viot,viol)
c         l2=lamo(2,tn,vrn,vtn,vln,tim,vior,viot,viol)
c
          qs=qs+(l1*par(1,i,ig)+l2*par(2,i,ig))*pari(1,i,ig)
    6     continue
          m=1
    4     if(m.gt.mass(11)) go to 3
            r1=alt*par(6,i,ig)**2
            r1=(qs+r1)*dts
            r2=(2*par(6,i,ig)+pari(1,i,ig))*alt
            r2=1+r2*dts
      if(r2.eq.0)print 901,par(6,i,ig),i,ig,pari(1,i,ig)
  901 format('r2=0',10g12.4)
            par(6,i,ig)=(moli+r1)/r2
            m=m+1
            go to 4
    3     continue
    2   continue
    1 continue
      return
      end

