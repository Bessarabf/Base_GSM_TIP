c ver 20/03/2019
c 4d massive dissociation rates (1 - 1/cm3/s; 2 - erg/cm3/s )
      subroutine r_dis(qdis,ano2,tem,gkoor,g,rads,solu,nsu,del,
     *           nh,its,ids,uts)

      USE mo_bas_gsm, ONLY: pi,om,bk,re,amO2
      dimension ano2(its,ids,nh),tem(its,ids,nh)
     *         ,g(nh),rads(nh)
     *         ,qdis(2,its,ids,nh),gkoor(2,its,ids)
 
      dimension solu(nsu),sp(12)
      data 
cc      cross sections F10.7= 70
cc   *     sp/0.12e-18,1.39e-18,13.2e-18,10.5e-18,2.39e-18,
cc   *        0.87e-18,0.28e-18,0.01e-18,0.,0.,0.,0./
c       cross sections F10.7=115
     *     sp/7.29E-19,4.30E-18,4.03E-19,4.69E-19,2.29E-18,9.40E-18,
     *        1.37E-17,1.01E-17,5.98E-18,2.55E-18,1.08E-18,3.93E-19/
c       cross sections F10.7=180
cc   *     sp/8.14e-19,4.27e-18,4.51e-19,4.69e-19,2.31e-18,9.48e-18,
cc   *        1.37e-17,1.01e-17,6.01e-18,2.58e-18,1.09e-18,3.92e-19/

      cr=pi/180.
      nsu05=nsu/2
!!!	qdis=1.e-10 !!!!!!!!!!!!!!!!
      sum=0.
      do i=2,its-1
         do j=1,ids

	     gshir=gkoor(1,i,j)*cr
           gdol=gkoor(2,i,j)*cr
           gshir=0.5*pi-gshir
c     !!!!!!!  zenith angle   !!!!!
           coshi=sin(gshir)*sin(del)+cos(gshir)*cos(del)*
     *           cos(om*(uts-43200.)+gdol)
           hi=acos(coshi)
           sumI=0.
           do k=nh,1,-1
              ra=sqrt(rads(k)*(rads(k)+2.*re))
              alfa=atan(re/ra)
              him=pi-alfa
              if(hi.le.him) then
      
                hO2=(bk*tem(i,j,k))/(amO2*g(k))
                reh=(rads(k)+re)/hO2
              !  Chepmen function
                chep=chept(reh,hi)
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                IF(K.EQ.NH) THEN 
                  sumI=anO2(i,j,k)*hO2
                else
                  sumI=sumI+0.5*(anO2(i,j,k)+anO2(i,j,k+1))*
     *                 (rads(k+1)-rads(k))    
                ! print sumI 
                end if
                sumL=0. 
                sumErg=0.                                               
                do l=1,nsu05                                   
                   tau=sp(l)*sumI*chep 
                   expTAU=exp(-tau)                                  
                   sumL=sumL+sp(l)*solu(l)*expTAU 
                   sumErg=sumErg+sp(l)*solu(l+nsu05)*expTAU                
                    ! print*,tau,k,l  
                end do                                  
                qdis(1,i,j,k)=sumL*1.e9 ! *anO2(i,j,k)!  
                qdis(2,i,j,k)=sumErg *anO2(i,j,k) 
              end if
            end do                           
          end do
      end do
      
      return 
      end

