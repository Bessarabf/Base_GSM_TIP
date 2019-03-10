      subroutine ionizu(sole,solen,rads,par,gkoor,ut,dolm,
     *       ddolgs,del,nh,its,ids,nse,kpars)
      dimension am(6),rads(nh),ai(4),sa(4,15),si(4,15),
     *       sole(nse),solen(nse),gkoor(2,its,ids),
     *       par(kpars,nh,its)
c     double precision tau0,xp,xxp,sx,reh,chep,derf,hsm,
c    *       re,pi,xlg,cr,xig,xi,cx,chab,ab,ex1,ex2
      data am/32.,28.,16.,30.,4.,1./,bk/1.38e-16/,
     *       ae/1.66e-24/,g0/980.665/,pi/3.1415926d0/,
     *       re/6371.02d5/,om/7.2722e-5/
      data si/4.63, 4.98, 2.49, 0.0,  7.0, 4.23, 6.0, 0.0 ,12.68, 9.4,
     *        7.20, 0.0 ,13.0 , 9.0,  7.5, 0.0 ,12.5, 9.36, 7.58, 0.0,
     *        13.0, 8.56, 9.12, 0.0, 15.5, 8.5 , 9.2, 0.0 ,19.53,18.00,
     *        10.5, 0.0 ,24.12,23.47,12.69,0.0 ,16.27,20.58,
     *        8.16, 0.0 , 5.93, 0.0 , 4.0 ,0.0 , 4.8,3*0.0 ,2.5,0.0,
     *        2*0.0,1.,0.0,5*0.,2.02/
      data sa/0.46, 0.34, 0.16, 0. , 1.67, 0.8, 1.00, 0. , 4.68 ,3.48,
     *        4.24, 0.  , 6.5 , 5.0, 5.0 , 0. , 8.5 , 5.4, 6.50, 0.0 ,
     *        12.0, 7.20, 8.00, 0.0, 15.5, 8.5, 9.2 , 0. ,19.00,18.50,
     *        10.5, 0.  ,24. , 24.4, 12.7, 0. ,27.3 ,27.3, 8.50, 0.0 ,
     *        16.2, 9.0 ,  4.0, 0. , 7.7 , 0.02, 0. , 0.0, 4.00, 0.0 ,
     *         0. , 0.  ,  1.6, 0.0, 0.0 , 0.0 , 0.01,0.0, 0.0 , 2.42/

      cr=180./pi
      bkg=bk/g0/ae
      k=4
      kl=nse
      nd=dolm/ddolgs+1
c
      do 1 ig=1,its
        rlat=gkoor(1,ig,nd)/cr
        rlat=pi/2.-rlat
        d olgg=gkoor(2,ig,nd)
        hl=ut+dolgg*3600./15.
        cx=sin(del) *sin(rlat)+cos(del)*cos(rlat)*cos(om*(hl-43200.))
        sx=sqrt(1.0-cx*cx)
        xi=acos(cx)
c       if(cx.lt.0.) xi=pi-xi
        xig=xi*cr
c ***** for Layman night radiation only
        if (xig.gt.90.0) then
         xlg=-0.96*cos(1.2*(xig-180.0)/cr)
        else
         xlg=0.0
        end if
c ***************************************
        do 9 i=1,k
          ai(i)=0.
    9   continue
        do 2 iv=1,nh
          i=nh-iv+1
        cone=0.
        sm=0.
          do 3 j=1,k
            cone=cone+par(j,i,ig)
            sm=sm+am(j)*par(j,i,ig)
    3     continue
          sm=sm/cone
          hsm=sm/(bkg*par(7,i,ig))
          reh=(rads(i)+re)*hsm
c ********* Chepmen function**********
          if(xig.LT.80.0) THEN
        	  chep=1.0/cx
          ELSE
	          chep=chept(reh,xi)
          END IF
          if(i.eq.nh) go to 7
            do 8 m=1,k
              ai(m)=ai(m)+(par(m,i,ig)+par(m,i+1,ig))*
     *        (rads(i+1)-rads(i))/2.
    8       continue
    7     continue
          do 4 j=1,k
            fs=0.
            do 5 l=1,kl
              ab=0.
              do 6 m=1,k
                ab=ab+sa(m,l)*ai(m)
    6         continue
              ab=ab*1.e-18
              chab=chep*ab
c             if(chab.ge.100.0)chab=100.0
              ex1=exp(-chab)
c             if(ab.ge.100.0)ab=100.0
              ex2=exp(-ab+xlg)
              if(ex1.le.1.e-30) ex1=1.e-30
              if(ex2.le.1.e-30) ex2=1.e-30
c
              fs=fs+si(j,l)*(sole(l)*ex1+solen(l)*ex2)
    5       continue
            fs=fs*1.e-9
            mm=j+12
            if(j.eq.3) mm=16
            if(j.eq.4) mm=15
            par(mm,i,ig)=par( j,i,ig)*fs
            if(par(mm,i,ig).le.1.e-20)  par(mm,i,ig)=1.e-20
    4     continue
    2   continue
               do 11 ii = 13 , 16
              a1=alog(par(ii,nh-1,ig))
              a2=alog(par(ii,nh-2,ig))
              a3=alog(par(ii,nh-3,ig))
              ac=a3+3.*a1-3.*a2
              par(ii,nh,ig)=exp(ac)
   11      continue
    1 continue
      return
      end

