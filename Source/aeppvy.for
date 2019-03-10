      subroutine aeppVY(mass,nl,idt,ddolgs,dtet,AL,Dst,ut,del0,E0,Q)
      !      nl=its, idt=ids
      dimension E0(nl,idt),Q(nl,idt),MASS(30)
      allocatable qe0(:,:),ebpp1(:),ebpp2(:),ebpp3(:),ebpp4(:),
     *            ebpp5(:),gmltm(:)
      allocate(qe0(nl,idt+1),ebpp1(idt+1),ebpp2(idt+1),ebpp3(idt+1),
     *         ebpp4(idt+1),ebpp5(idt+1),gmltm(idt+1))      
      data pi/3.14159265359/
      dpi=pi+pi
      if(AL.ge.-5.)AL=-5.
      do l=1,idt+1
	 
        l00=l
	  if(l.eq.idt+1)l00=1 

        dolg=(l-1)*ddolgs
        fi=dolg/180.*pi
        call magsm(ut,del0,fi,phism,1)
        tau=pi+phism
        if(tau.ge.dpi)tau=tau-dpi
        if(tau.lt.0.)tau=tau+dpi
        gmlt=tau/pi*12.
        gmltm(l)=gmlt
          if(gmlt.eq.0.)then
            b1e=64.64+0.01*al+1.25*1.e-6*al*al+.02*dst
            b2e=66.69+0.009*al+7.78*1.e-7*al*al+.022*dst
            b5e=70.89-0.002269047619*al+0.00175*dst
            b6e=71.74-0.000994*al-0.0038*dst
            b2i=66.25158347+0.005469936384*al-5.430946029e-7*al*al
     *      +0.02255*dst+0.361
            b4s=67.48+0.007*al+1.36e-6*al*al+0.048*dst
          end if
          if(gmlt.gt.0..and.gmlt.le.3.)then
            b1e=64.77+0.012*al+1.2*1.e-6*al*al+.03*dst
            b2e=67.1+0.01*al+.58*1.e-6*al*al+.025*dst
            b5e=73.2+0.003*al+.0026*dst
            b6e=74.2+0.004*al+.01*dst
            b2i=68.7886781+0.0060622409195*al-2.7154730145e-7*al*al
     *      +0.03035*dst
c            b4s=(71.25629735+0.006298636364*al-2.916666667e-7*al*al
c     *      -0.3755-0.00545*dst+67.48+0.007*al+1.36e-6*al
c     *      *al+0.048*dst)*.5
            b4s=70.19100384905+0.0064990462755*AL+5.3379164115E-007*AL
     *          *AL+0.024404761905*Dst
          end if
          if(gmlt.gt.3..and.gmlt.le.6.)then
            b1e=64.90363095+0.01525*al+6.547619048*1.e-6*al*al
     *      +0.056*dst
c            b2e=68.29+0.014*al+7.187*1.e-6*al*al+0.033625*dst
            b2e=67.69140422+0.01406017316*AL+7.187229437e-006*AL*AL
     *          +0.03076785714*Dst+0.1494642857
            b5e=75.85+0.0061*al+0.0043*dst
            b6e=77.2+0.00746*al+0.012*dst
            b2i=70.06127273+0.006654545455*al+0.03815*dst+0.9035
            b4s=71.25629735+0.006298636364*al-2.916666667e-7*al*al
     *          +0.03125595238*Dst+0.7564880952
c     *      -0.3755-0.00545*dst
          end if
          if(gmlt.gt.6..and.gmlt.le.9.)then
            ecps=66.129+0.0142*al+6.871*1.e-6*al*al+.0524*dst
            pcps=72.6485+0.01603396486*al+8.97319837*1.e-6*al*al+
     *      0.023827*dst
            ebps=73.35+0.01145*al+3.292*1.e-6*al*al+0.02354761905*dst
            pbps=77.71435386+0.01015437229*al+2.79004329*1.e-6*al*al
     *      +0.01985714286*dst
            pllb=78.99303084381+0.01654489177*al+5.851731602*1.e-6*al
     *      *al+0.02108571429*dst
          end if
          if(gmlt.gt.9..and.gmlt.le.12.)then
            ecps=67.21+0.006470909091*al+0.0219*dst 
            pcps=74.03+0.00579197*al-3.40909*1.e-8*al*al+0.0293*dst 
            ebps=76.2+0.00806*al+0.0295*dst 
            pbps=78.11+0.007922727273*al+0.0186*dst
            pllb=80.53+0.01210984848*al+2.276515152*1.e-6*al*al+0.0293
     *      *dst
          end if
          if(gmlt.gt.12..and.gmlt.le.15.)then
            ecps=69.81670238+0.008624285714*al+2.071428571*1.e-6*al*al
     *      +0.04185*dst 
            pcps=74.2+0.01015458874*al+2.336580087*1.e-6*al*al
     *      +0.04185*dst 
            ebps=77.97+0.0137058807*al+3.324815877*1.e-6*al*al
     *      +0.0511*dst   
            pbps=79.984+0.01357380952*al+4.130952381*1.e-6*al*al
     *      +0.04585*dst
            pllb=81.18+0.01602409091*al+4.526515152*1.e-6*al*al
     *      +0.04745*dst
          end if
          if(gmlt.gt.15..and.gmlt.le.18.)then
            ecps=70.22+0.01153809524*al+5.345238095*1.e-6*al*al
     *      +0.06485*dst 
            pcps=73.266+0.01364214876*al+6.432113341*1.e-6*al*al
     *      +0.04532857143*dst 
            ebps=74.4+0.01529090909*al+6.735930736*1.e-6*al*al
     *      +0.0407*dst   
            pbps=78.+0.01693571429*al+8.94047619*1.e-6*al*al
     *      +0.02945*dst
            pllb=80+0.02479285714*al+1.272619048*1.e-5*al*al
     *      +0.02945*dst
          end if
          if(gmlt.gt.18..and.gmlt.le.21.)then
            b1e=68.91+0.012*al+3.85*1.e-6*al*al+0.05691428571*dst
c            b2e=71.14+0.011*al+3.61*1.e-6*al*al+0.05691428571*dst
            b2e=69.73307168+0.01133142191*AL+3.60955711E-006*AL*AL
     *          +0.9166517857+0.04044642857*Dst
            b5e=74.93+0.01275247934*al+8.820936639*1.e-6*al*al
     *      +0.0229*dst
            b6e=76.52+0.01780134525*al+1.155390033*1.e-5*al*al
     *      +0.01865*dst
            b2i=68.22612587+0.009239067599*al+2.610722611e-6*al*al
     *      +0.04798571429*dst+1.174142857
            b4s=71.48292992+0.01292863636*al+5.65530303e-6*al*al
     *      +0.024*dst+0.3
          end if
          if(gmlt.gt.21..and.gmlt.le.24.)then
            b1e=64.64+0.01*al+1.25*1.e-6*al*al+.02*dst
c            b2e=66.69+0.009*al+7.78*1.e-7*al*al+.022*dst
            b2e=66.71830828+0.009183808742*AL+7.78062795e-7*AL*AL
     *          +0.0226*Dst+0.07266666667
            b5e=70.89-0.002269047619*al+0.00175*dst
            b6e=71.74-0.000994*al-0.0038*dst
            b2i=66.25158347+0.005469936384*al-5.430946029e-7*al*al
     *      +0.02255*dst+0.361
c            b4s=67.48+0.007*al+1.36e-6*al*al+0.048*dst
            b4s=68.39+0.007*AL+1.36e-6*AL*AL+0.017*Dst
          end if
          if(gmlt.ge.0..and.gmlt.le.6..or.gmlt.gt.18..
     *    and.gmlt.le.24.)then
            ebpp1(l)=b1e
            ebpp2(l)=b2e
            ebpp3(l)=b4s
            ebpp4(l)=b5e
            ebpp5(l)=b6e
          end if
          if(gmlt.gt.6..and.gmlt.le.18.)then
            ebpp1(l)=ecps
            ebpp2(l)=pcps
            ebpp3(l)=ebps
            ebpp4(l)=pbps
            ebpp5(l)=pllb
          end if
        call aurprecip(gmlt,al,f1,f2,
     *         f4,f5,e1,e2,e4,e5,eaop,
     *         edaz,esdp,faop,fdaz,fsdp)
        do j=1,nl/2+1
          tet=(j-1)*dtet

          gmlat=90.-tet
          if(gmlt.ge.0..and.gmlt.le.6..or.gmlt.gt.18..
     *       and.gmlt.le.24.)then
            if(gmlat.gt.b1e.and.gmlat.le.b2e)then
              qe0(j,l)=f1
              E0(j,l00)=e1
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.b2e.and.gmlat.le.b4s)then
              qe0(j,l)=f2
              E0(j,l00)=e2
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.b4s.and.gmlat.le.b5e)then
              qe0(j,l)=f4
              E0(j,l00)=e4
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.b5e.and.gmlat.le.b6e)then
              qe0(j,l)=f5
              E0(j,l00)=e5
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
          end if
          if(gmlt.gt.6..and.gmlt.le.18.)then
            if(gmlat.gt.ecps.and.gmlat.le.pcps)then
              qe0(j,l)=fdaz
              E0(j,l00)=edaz
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.ge.ebps.and.gmlat.le.pbps)then
              qe0(j,l)=faop
              E0(j,l00)=eaop
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
            if(gmlat.gt.pbps.and.gmlat.le.pllb)then
              qe0(j,l)=fsdp
              E0(j,l00)=esdp
              Q(j,l00)=qe0(j,l)/E0(j,l00)*6.24e8
            end if
          end if

        end do
      end do
! south hemisphere 
      
      do j=1,nl/2+1
        jc=nl-j+1
	  do l=1,idt+1
           l00=l
           if(l.eq.idt+1)l00=1 
           E0(jc,l00)=E0(j,l00)
           Q(jc,l00)=Q(j,l00)
        end do
      end do
		  
        if(mass(30).eq.1)then	
          open(1,file='epres',access='direct',recl=1371)
          nrec=nrec+1
          write(1,rec=nrec)ut
	      nrec=nrec+1
          write(1,rec=nrec)AL
	      nrec=nrec+1
          write(1,rec=nrec)Dst
	      do j=1,nl/2+1
	        nrec=nrec+1
	        write(1,rec=nrec)(E0(j,l),l=1,idt)
	      end do
	      do j=1,nl/2+1
	        nrec=nrec+1
	        write(1,rec=nrec)(qe0(j,l),l=1,idt)
	      end do
	      do j=1,nl/2+1
	        nrec=nrec+1
	        write(1,rec=nrec)(Q(j,l),l=1,idt)
	      end do
          close(1)
          open(1,file='bound',access='direct',recl=147)
	      nrec=nrec+1
          write(1,rec=nrec)ut
	      nrec=nrec+1
          write(1,rec=nrec)AL
	      nrec=nrec+1
          write(1,rec=nrec)Dst
	      do l=1,idt
	        nrec=nrec+1
            write(1,rec=nrec)gmltm(l),ebpp1(l),ebpp2(l),ebpp3(l),
     *          ebpp4(l),ebpp5(l)
	      end do
          close(1)
        end if

      print*,' '
      print 100,(ebpp1(i),i=1,idt)
      print*,' '
  100 format(' ',5f9.3)
      deallocate(qe0,ebpp1,ebpp2,ebpp3,ebpp4,ebpp5,gmltm)      
      return
      end



      subroutine aurprecip(gmlt,al,f1,f2,
     *       f4,f5,e1,e2,e4,e5,eaop,edaz,esdp,
     *       faop,fdaz,fsdp)
      if(gmlt.eq.0.)then
        f1=exp(0.6571873362*alog(abs(al))-3.393664639) 
        f2=exp(0.709435388*alog(abs(al))-2.682839591)
        f4=exp(0.4970209835*alog(abs(al))-0.8709023689)
        f5=exp(0.1248129607*alog(abs(al))-1.84467514)
        e1=exp(0.2783876967*alog(abs(al))-0.553631512) 
        e2=exp(0.3491325617*alog(abs(al))-0.5661737632)
        e4=exp(0.2860529128*alog(abs(al))-0.6452854749)
        e5=exp(0.377738*alog(abs(al))-2.183403144)
      end if
      if(gmlt.gt.0..and.gmlt.le.3.)then
        f1=exp(0.5010929633*alog(abs(al))-2.057026824)
        f2=exp(0.7668563524*alog(abs(al))-2.961112584)
        f4=exp(0.5390297179*alog(abs(al))-1.599382872)
        f5=4.254511278E-5*abs(AL)+0.2174526316
c        f5=abs(f5) !!! correction 28.08.2015
        faop=exp(0.7229942079*alog(abs(AL))-2.72199829)
        e1=exp(0.208705056*alog(abs(al))+0.1701497954)
        e2=exp(0.2827463961*alog(abs(al))-0.2616186626)
        e4=exp(0.2847802859*alog(abs(al))-1.017263627)
        e5=exp(0.2500022797*alog(abs(al))-1.712236379)
        eaop=0.8683157161*alog(abs(AL))-2.02028153
      end if
      if(gmlt.gt.3..and.gmlt.le.6.)then
c        f1=exp(0.5010929633*alog(abs(al))-2.057026824)
c        f2=exp(0.7668563524*alog(abs(al))-2.961112584)
c        f4=exp(0.5390297179*alog(abs(al))-1.599382872)
c        f5=4.254511278e-5*al+0.2174526316
c        f5=abs(f5) !!! correction 28.08.2015
        f2=exp(0.8159742601*alog(abs(AL))-3.195652962)
        f4=exp(0.7038019846*alog(abs(AL))-3.54303041)
        f1=exp(0.4158926226*alog(abs(AL))-1.268071995)
        f5=exp(0.1468195661*alog(abs(AL))-2.846272624)
        faop=exp(0.8862633552*alog(abs(AL))-3.878830031)
        e1=exp(0.1667090503*alog(abs(al))+0.6541825456)
        e2=exp(0.198093427*alog(abs(al))+0.1380986798)
        e4=exp(0.2812593182*alog(abs(al))-1.615408825)
        e5=exp(-0.0523166127*alog(abs(al))-0.3659097527)
        eaop=0.7665735864*alog(abs(AL))-1.723066197
      end if
      if(gmlt.gt.6..and.gmlt.le.9.)then
        eaop=exp(0.08592111476*alog(abs(al))-0.2500901102)
        edaz=exp(0.1915357972*alog(abs(al))+0.5662103641)
        esdp=-5.e-5*abs(al)+0.3158333333
        faop=exp(0.4138865686*alog(abs(al))-1.567065282)
        fdaz=exp(0.3313262308*alog(abs(al))-1.443740347)
        fsdp=0.0007166666667*abs(al)+0.4697222222
c        fsdp=abs(fsdp) !!! correction 28.08.2015
      end if
      if(gmlt.gt.9..and.gmlt.le.12.)then
        eaop=exp(0.08496450811*alog(abs(al))-0.5313937638) 
        edaz=exp(0.2007310553*alog(abs(al))+0.7161005602)
        esdp=-1.e-5*abs(AL)+0.2756111111
        faop=exp(0.1615595083*alog(abs(al))-0.7343816142)
        fdaz=exp(0.07503442233*alog(abs(al))-0.7564862408)
c        fsdp=0.0006016666667*al+0.5025833333
c        fsdp=abs(fsdp) !!! correction 28.08.2015
cc       fsdp0=0.2638098631*alog(abs(al))-1.826307057
        fsdp=0.1815317224*alog(abs(AL))-0.2868825239
      end if
      if(gmlt.gt.12..and.gmlt.le.15.)then
        eaop=exp(0.06244459649*alog(abs(al))-0.8492598268) 
        edaz=exp(0.04497766381*alog(abs(al))+1.217024722)  
        esdp=5.142857143e-5*abs(al)+0.2312380952
        faop=exp(0.3052293503*alog(abs(al))-1.577706858)  
        fdaz=exp(-0.145536954*alog(abs(al))-1.148174158)
c        fsdp=0.00082*al+0.409
c        fsdp=abs(fsdp) !!! correction 28.08.2015
        fsdp=exp(0.3548040564*alog(abs(AL))-2.341480608)
      end if
      if(gmlt.gt.15..and.gmlt.le.18.)then
        eaop=exp(0.2654109679*alog(abs(al))-1.351216348) 
        edaz=exp(-0.05754892011*alog(abs(al))+1.149209589)  
        esdp=-2.285714286e-5*abs(al)+0.3101904762
        faop=exp(0.5976870883*alog(abs(al))-2.775956307)  
        fdaz=exp(-0.1464972011*alog(abs(al))-1.224727104) 
        fsdp=0.001574285714*abs(al)+0.3127142857
c        fsdp=abs(fsdp)  !!! correction 28.08.2015
      end if
      if(gmlt.gt.18..and.gmlt.le.21.)then
        f1=exp(0.5024018313*alog(abs(al))-3.406987641)
        f2=exp(0.5576887485*alog(abs(al))-3.099012084)
        f4=exp(0.7512055891*alog(abs(al))-2.799803464)
        f5=exp(0.4294515713*alog(abs(al))- 3.602047784)
        faop=exp(0.806511179*alog(abs(AL))-3.313109732)
        e1=exp(0.1264477554*alog(abs(al))-0.265670415)
        e2=exp(0.1617507433*alog(abs(al))-0.28674081)
        e4=exp(0.4355726882*alog(abs(al))-1.732201529)  
        e5=exp(0.5313973653*alog(abs(al))-3.340605919)
        eaop=0.385046894*alog(abs(AL))-1.365733462
      end if
      if(gmlt.gt.21..and.gmlt.le.24.)then
        f1=exp(0.6571873362*alog(abs(al))-3.393664639) 
        f2=exp(0.709435388*alog(abs(al))-2.682839591)
        f4=exp(0.4970209835*alog(abs(al))-0.8709023689)
        f5=exp(0.1248129607*alog(abs(al))-1.84467514)
        faop=3.504970748*alog(abs(AL))-13.14906522
        e1=exp(0.2783876967*alog(abs(al))-0.553631512) 
        e2=exp(0.3491325617*alog(abs(al))-0.5661737632)
        e4=exp(0.2860529128*alog(abs(al))-0.6452854749)
        e5=exp(0.377738*alog(abs(al))-2.183403144)
        eaop=1.029401751*alog(abs(AL))-2.632421127
      end if
      return
      end   