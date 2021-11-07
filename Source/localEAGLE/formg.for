      subroutine formg(i,gg,gk,z,om,td,p,wt2,wt3,swt,cwt,df,dfs,ap4)
      dimension gk(45,6),g(45),p(12)
      do 4 j=1,45
    4 g(j)=gk(j,i)
      gk12=g(12)
      da=ap4+(g(13)-1.)*(ap4+(exp(-gk12*ap4)-1.)/gk12)
      if(i.le.3.) go to 3
         da=da*g(42)
         if(z.lt.140.) go to 1
            zz=z
            go to 2
    1       zz=140.
    2    continue
         arg1=exp(g(45)*(150.-zz))
         cwt=cwt*(1.+g(43)*arg1)
         swt=swt*(1.+g(44)*arg1)
    3 arg1=cos(om*(td-g(21)))
      arg2=1.+g(9)*dfs
      om2=om+om
      gg=1.+g(3)*p(2)+g(4)*p(4)+g(2)*p(1)+g(5)*df+g(6)*df*df+
     *g(7)*dfs+(g(10)+g(11)*p(2))*da
      gg=gg+g(14)*cos(om*(td-g(15)))
      gg=gg+(g(16)+g(17)*p(2))*cos(om2*(td-g(18)))
      gg=gg+(g(19)*p(1)+g(20)*p(3))*(1.+g(8)*dfs)*arg1
      gg=gg+g(22)*p(1)*cos(om2*(td-g(23)))
      gg=gg+(g(24)*p(5)+g(25)*p(7)+g(26)*p(8)+(g(27)*p(5)+g(28)*
     *p(6))*arg1)*arg2*cwt
      gg=gg+(g(29)*p(5)+g(30)*p(7)+g(31)*p(8)+(g(32)*p(5)+g(33)*
     *p(6))*arg1)*arg2*swt
      gg=gg+(g(34)*p(9)+g(35)*p(11)+g(36)*p(10)*arg1)*arg2*cos(wt2)
      gg=gg+(g(37)*p(9)+g(38)*p(11)+g(39)*p(10)*arg1)*arg2*sin(wt2)
      gg=gg+p(12)*arg2*(g(40)*cos(wt3)+g(41)*sin(wt3))
      return
      end

