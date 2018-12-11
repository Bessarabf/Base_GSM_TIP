c         p(1)=p10      p(7)=p31
c         p(2)=p20      p(8)=p51
c         p(3)=p30      p(9)=p22
c         p(4)=p40      p(10)=p32
c         p(5)=p11      p(11)=p42
c         p(6)=p21      p(12)=p33
c
      subroutine polinm(tet,p)
      dimension p(12)
      c=cos(tet)
      s=sin(tet)
      c2=c*c
      s2=s*s
      p(1)=c
      p(2)=(3.*c2-1.)*0.5
      p(3)=(5.*c2-3.)*c*0.5
      p(4)=((7.*c2-6.)*5.*c2+3.)*0.125
      p(5)=s
      p(6)=3.*s*c
      p(7)=(5.*c2-1.)*s*1.5
      p(8)=((3.*c2-2.)*7.*c2+1.)*s*1.875
      p(9)=3.*s2
      p(10)=15.*s2*c
      p(11)=(7.*c2-1.)*s2*7.5
      p(12)=15.*s2*s
      return
      end
 
