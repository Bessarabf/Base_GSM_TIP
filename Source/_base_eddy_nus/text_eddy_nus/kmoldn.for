      subroutine k moldn(dp,an1,an2,an3,an6,ro)
c
c     . . p/p ras. k. mol. dif. to trexkomponent.
c     sredi   v odnoj tocke
      dimension
     *          dp(3,3)
      data am1,am2,am3 /53.12e-24,46.51e-24,26.56e-24/,
     *     bk,pi/1.38e-16,3.1415296/,
     *     d1,d2,d3/3.6e-8,3.7e-8,2.e-8/
      coef=3./8.*sqrt(0.5*bk/pi)
      sum=an1+an2+an3
c     . . . otnosheie mol mass
      otm12=am1/am2
      otm13=am1/am3
      otm23=am2/am3
c     . . . srednij diametr vzaim.
      sig12=0.5*(d1+d2)
      sig13=0.5*(d1+d3)
      sig23=0.5*(d2+d3)
      sig12=sig12*sig12
      sig13=sig13*sig13
      sig23=sig23*sig23
c     . . . obrat. massa
      obr1=1./am1
      obr2=1./am2
      obr3=1./am3
      b12=sqrt(obr1+obr2)
      b13=sqrt(obr1+obr3)
      b23=sqrt(obr2+obr3)
      sqt=sqrt(an6)
      value=coef*sqt
      d12=value/sig12*(b12/sum)
      d13=value/sig13*(b13/sum)
      d23=value/sig23*(b23/sum)
      d32=d23
      ro1=an1*am1
      ro2=an2*am2
      ro3=an3*am3
      znam=(an1*d23+an2*d13+an3*d12)*ro
      rez=sum/znam
      dp(1,1)=d13*d12*(ro2+ro3)*rez
      dp(1,2)=-d23*d12*ro1*rez/otm12
      dp(1,3)=-d23*d13*ro1*rez/otm13
      dp(2,1)=-d12*d13*ro2*rez*otm12
      dp(2,2)=d23*d12*(ro1+ro3)*rez
      dp(2,3)=-d23*d13*ro2*rez/otm23
      dp(3,1)=-d12*d13*ro3*rez*otm13
      dp(3,2)=-d23*d12*ro3*rez*otm23
      dp(3,3)=d23*d13*(ro2+ro1)*rez
      return
      end
