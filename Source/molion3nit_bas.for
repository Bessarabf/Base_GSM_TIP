!!!!  ver 11.03.2019 - O2+ and NO+ write in file4
!!!!  NEWTON ITERATION
      subroutine molio3Nit_bas(par,ids,its,nh,kpars,pari,ins,
     *                  mass,dts,j,ntr,key) 
	                 
      dimension par(kpars,nh,its),
     *          pari(ins,nh,its),mass(30),al(10)

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!TABLE. The list of chemical reactions.
!1.		$O_{2}^{+}+e\to O+O$
!2.		$O_2^++NO \to NO^++O_2$
!3.		$O_{2}^{+}+N\to N{{O}^{+}}+O$
!4.		$O_{2}^{+}+{{N}_{2}}\to N{{O}^{+}}+NO$
!5.		$N{{O}^{+}}+e\to N+O$
!6.		${{O}^{+}}+{{O}_{2}}\to O_{2}^{+}+O$
!7.		${{O}^{+}}+{{N}_{2}}\to N{{O}^{+}}+N$
!8.		$N_{2}^{+}+O\to N{{O}^{+}}+N$
!9.		$N_{2}^{+}+{{O}_{2}}\to O_{2}^{+}+{{N}_{2}}$
!10.		$N_{2}^{+}+e\to N+N$

! reaction rates (first approcsimation)     
      
      data al/2.25e-7,6.3E-10,1.8E-10,1.0E-15,4.5e-7,
     *        2.0E-11,1.2E-12,1.4E-10,6.0E-11,4.3e-7/	
! 
! cO2I = O2+ ions on previos TIME level
! cNOI = NO+ ions on previos TIME level
! cO2pl(nh,its,ids) = massiv O2+ and 
! cNOpl(nh,its,ids) = NO 
! cn1I = O2+
! cn2I = NO+
! cn3I = O+
! cn4I = N2+
! 
! cn1I,cn2I etc - previos iteration
! cn1I_N,... - Next iteration  
! dn1I, dn2I - delta for O2+ & NO+ concentrations 
!


	do ig=1,its-1 
	  do k=1,nh
!!!!!!!!! initialisation 
            cO2=par(1,k,ig)    ! cO2	     
            cN2=par(2,k,ig)    ! cNO 
            cO=par(3,k,ig)     ! cO 
	    cNO=par(4,k,ig)    ! cNO
            if(mass(21).eq.0) then 
	      cN=0.
	    else
              cN=par(5,k,ig)
            end if
	    te=par(9,k,ig)
    ! ionization rates
	    qO2=par(13,k,ig) 
	    qN2=par(14,k,ig)
	    qNO=par(15,k,ig)
	    qO=par(16,k,ig)  
            Q3= qO  
    ! lost 
            pL3=al(6)*cO2+al(7)*cN2
    !     O+       !!!!!!!! 
           if (k.lt.ntr) then
                  cn3I=Q3/pL3     
           else 
                  cn3I=pari(1,k,ig)
           end if
           cNe=par(6,k,ig)  ! Ne значения на i-ом временном шаге 
          ! first step before iteration O2+ & NO+
 !             key=0
           if (key.eq.0) then
              cn1I=.5*par(6,k,ig)! Q1(k)/pL1(k)
              cn2I=.5*par(6,k,ig)! Q2(k)/pL2(k)
              par(18,k,ig)=cn2I
              par(19,k,ig)=cn2I
           else
              cn1I=par(18,k,ig) !!cO2pl(k,ig,j)
              cn2I=par(19,k,ig) !!!cNOpl(k,ig,j)
	   end if
                ! lost rates
           pL4=al(8)*cO+al(9)+cO2+al(10)*cNe
                ! source 
           Q4= qN2
           cn4I=Q4/pL4
           Q1= qO2+al(6)*cn3I*cO2+al(9)*cn4I*cO2
             ! Newton iteration
           dnorm=11.
           it=0
           do while(dnorm.gt.0.001)
	!     do it=1,3 
              pL1=al(1)*cNe+al(2)*cNO+al(3)*cN+al(4)*cN2
              pL2=al(5)*cNe
              Q2=qNO+(al(2)*cNO+al(3)*cN+al(4)*cN2)*cn1I+
     *           al(7)*cn3I*cN2+al(8)*cn4I*cO
              F1=cn1i*(1+pL1*dts)-par(18,k,ig)-Q1*dts
              F2=cn2i*(1+pL2*dts)-par(19,k,ig)-Q2*dts
              ! Jacobian  
              G11=1.+pL1*dts+al(1)*cn1I*dts
              G12=0.
              G21=(al(2)*cNO+al(3)*cN+al(4)*cN2)*dts
              G22=1.+pL2*dts+al(5)*cn2I*dts
              ! delta O2+ & NO+
              dn1I=-F1/G11
              dn2I=-(F2+G21*F1/G11)/G22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !dnorm=sqrt(dn1I**2+dn2I**2)
              dnorm=abs(dn1I)/cn1I
              dnorm2=abs(dn2I)/cn2I
              if (dnorm2.gt.dnorm) dnorm=dnorm2
              cn1I=cn1I+dn1I
              cn2I=cn2I+dn2I
              cNe=cn1I+cn2I+cn3I+cn4I
              it=it+1
              if(it.ge.6) EXIT
           end do
!!!!   correction negative values
           if(cn1I.lt.0.) cn1I=par(18,k-1,ig)
           if(cn2I.lt.0.) cn2I=par(19,k-1,ig)   
           par(18,k,ig)=cn1I
           par(19,k,ig)=cn2I
 !          write(10,*) ig,k,j,cO2pl(k,ig,j),cNOpl(k,ig,j),dnorm,it
           par(6,k,ig)=cn1I+cn2I+cn4I !!!!cO2pl(k,ig,j)+cNOpl(k,ig,j)
         end do           
	end do
	return
	end
        
        