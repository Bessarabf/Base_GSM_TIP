	subroutine aepp(mass,Kp,nl,idt,ddolgs,dtet,ut,del0,E0,Q,its,ids)
     
      dimension Af(6),Bf(6),Cf(6),Df(6),Ae(6),Be(6),Ce(6),De(6)

      dimension E0(its,ids),Q(its,ids),mass(30)
     	
      dimension Acf(6,6),Bcf(6,6),Ccf(6,6),Dcf(6,6),Asf(6,6),Bsf(6,6),
     *Csf(6,6),Dsf(6,6),Ace(6,6),Bce(6,6),Cce(6,6),Dce(6,6),Ase(6,6),
     *Bse(6,6),Cse(6,6),Dse(6,6)
      dimension Akf(6),Bkf(6),Ckf(6),Dkf(6),Ake(6),Bke(6),Cke(6),Dke(6)
      real Kp,Kpm1,Kpm2,Kp_model(6)

	allocatable  qe0(:,:)
     	allocate (qe0(its,ids))

	data pi/3.14159265359/
      data Kp_model/.75,2.25,3.75,5.25,7.,9./
	
      open(1,file='Zhang_Paxton.dat')
      read(1,1)(Akf(i),i=1,6)
      do j=1,6
        read(1,1)(Acf(j,i),i=1,6)
        read(1,1)(Asf(j,i),i=1,6)
      end do
      read(1,1)(Bkf(i),i=1,6)
      do j=1,6
        read(1,1)(Bcf(j,i),i=1,6)
        read(1,1)(Bsf(j,i),i=1,6)
      end do
      read(1,1)(Ckf(i),i=1,6)
      do j=1,6
        read(1,1)(Ccf(j,i),i=1,6)
        read(1,1)(Csf(j,i),i=1,6)
      end do
      read(1,1)(Dkf(i),i=1,6)
      do j=1,6
        read(1,1)(Dcf(j,i),i=1,6)
        read(1,1)(Dsf(j,i),i=1,6)
      end do
      read(1,1)(Ake(i),i=1,6)
      do j=1,6
        read(1,1)(Ace(j,i),i=1,6)
        read(1,1)(Ase(j,i),i=1,6)
      end do
      read(1,1)(Bke(i),i=1,6)
      do j=1,6
        read(1,1)(Bce(j,i),i=1,6)
        read(1,1)(Bse(j,i),i=1,6)
      end do
      read(1,1)(Cke(i),i=1,6)
      do j=1,6
        read(1,1)(Cce(j,i),i=1,6)
        read(1,1)(Cse(j,i),i=1,6)
      end do
      read(1,1)(Dke(i),i=1,6)
      do j=1,6
        read(1,1)(Dce(j,i),i=1,6)
        read(1,1)(Dse(j,i),i=1,6)
      end do
      close(1)
    1 format(6f12.7)
      dpi=pi+pi
c     Определение интервала для Kp индекса и Kpm1 и Kpm2
      if(Kp.le.0.75)then
        Kpm1=0.75
        Kpm2=2.25
	end if
      if(Kp.gt.0.75.and.Kp.le.9.)then
	  do i=1,5
	    if(Kp_model(i).lt.Kp.and.Kp_model(i+1).ge.Kp)then
		  Kpm1=Kp_model(i)
		  Kpm2=Kp_model(i+1)
	    end if
	  end do
	end if
	if(Kp.gt.9.)then
	  Kpm1=7.00
	  Kpm2=9.00
	end if
c     Нахождение Hp индекса и Hpm1 и Hpm2
	if(Kp.le.5.)Hp=38.66*exp(0.1967*Kp)-33.99
	if(Kp.gt.5.)Hp=4.592*exp(0.4731*Kp)+20.47
	if(Kpm1.le.5.)Hpm1=38.66*exp(0.1967*Kpm1)-33.99
	if(Kpm1.gt.5.)Hpm1=4.592*exp(0.4731*Kpm1)+20.47
	if(Kpm2.le.5.)Hpm2=38.66*exp(0.1967*Kpm2)-33.99
	if(Kpm2.gt.5.)Hpm2=4.592*exp(0.4731*Kpm2)+20.47
	if(Kpm1.eq.0.75)then
        i=1
        ip=2 
      end if
	if(Kpm1.eq.2.25)then
        i=2
        ip=3 
      end if
	if(Kpm1.eq.3.75)then
        i=3
        ip=4 
      end if
	if(Kpm1.eq.5.25)then
        i=4
        ip=5 
      end if
	if(Kpm1.eq.7.00)then
        i=5
        ip=6 
      end if
c     Определение углов
        do j=1,nl
          tet=(j-1)*dtet
          alat=90.-tet
          x=90.-abs(alat)  
          do l=1,idt
            if(j.gt.10.and.j.lt.28)then
              E0(j,l)=1.e-3
		  qe0(j,l)=0.
              Q(j,l)=0.  
              goto 2 
            end if  
            dolg=(l-1)*ddolgs
            fi=dolg/180.*pi
            call magsm(ut,del0,fi,phism,1)
            tau=pi+phism
            if(tau.ge.dpi)tau=tau-dpi
            if(tau.lt.0.)tau=tau+dpi
c     Нахождение коэффициентов Af1-Df1 Af2-Df2 Ae1-De1 Ae2-De2
            do k=i,ip
              Af(k)=Akf(k)
              Bf(k)=Bkf(k)
              Cf(k)=Ckf(k)
              Df(k)=Dkf(k) 
              Ae(k)=Ake(k)
              Be(k)=Bke(k)
              Ce(k)=Cke(k)
              De(k)=Dke(k)
              do m=1,6
                Af(k)=Af(k)+Acf(k,m)*cos(m*tau)+Asf(k,m)*sin(m*tau) 
                Bf(k)=Bf(k)+Bcf(k,m)*cos(m*tau)+Bsf(k,m)*sin(m*tau) 
                Cf(k)=Cf(k)+Ccf(k,m)*cos(m*tau)+Csf(k,m)*sin(m*tau) 
                Df(k)=Df(k)+Dcf(k,m)*cos(m*tau)+Dsf(k,m)*sin(m*tau) 
                Ae(k)=Ae(k)+Ace(k,m)*cos(m*tau)+Ase(k,m)*sin(m*tau) 
                Be(k)=Be(k)+Bce(k,m)*cos(m*tau)+Bse(k,m)*sin(m*tau) 
                Ce(k)=Ce(k)+Cce(k,m)*cos(m*tau)+Cse(k,m)*sin(m*tau) 
                De(k)=De(k)+Dce(k,m)*cos(m*tau)+Dse(k,m)*sin(m*tau)
              end do 
            end do
c     Нахождение Eom1 E0m2 Qm1 Qm2
           E0m1=Ae(i)*exp((x-Be(i))/Ce(i))/(1.+exp((x-Be(i))/De(i)))**2
            E0m2=Ae(ip)*exp((x-Be(ip))/Ce(ip))/(1.+exp((x-Be(ip))/
     *        De(ip)))**2
           Qm1=Af(i)*exp((x-Bf(i))/Cf(i))/(1.+exp((x-Bf(i))/Df(i)))**2
            Qm2=Af(ip)*exp((x-Bf(ip))/Cf(ip))/(1.+exp((x-Bf(ip))/
     *        Df(ip)))**2
c     Нахождение E0 - энергии
		f1=(Kpm2-Kp)/(Kpm2-Kpm1)
		f2=(Kp-Kpm1)/(Kpm2-Kpm1)
		E0(j,l)=f1*E0m1+f2*E0m2

c     Нахождение Q - потока
		f1=(Hpm2-Hp)/(Hpm2-Hpm1)
		f2=(Hp-Hpm1)/(Hpm2-Hpm1)
		qe0(j,l)=f1*Qm1+f2*Qm2
            Q(j,l)=qe0(j,l)/E0(j,l)*6.24e8 
    2     continue
	  end do
	end do
        if(mass(30).eq.1)then	
          open(1,file='epres',access='direct',recl=723)
	    nrec=nrec+1
          write(1,rec=nrec)ut
	    nrec=nrec+1
          write(1,rec=nrec)Kp
	    nrec=nrec+1
          write(1,rec=nrec)Hp
	    do j=1,10
	      nrec=nrec+1
	      write(1,rec=nrec)(E0(j,l),l=1,24)
	    end do
	    do j=1,10
	      nrec=nrec+1
	      write(1,rec=nrec)(qe0(j,l),l=1,24)
	    end do
	    do j=1,10
	      nrec=nrec+1
	      write(1,rec=nrec)(Q(j,l),l=1,24)
	    end do
        end if
      close(1)
	deallocate (qe0)
	return 
	end