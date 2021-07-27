      subroutine bonvec1(vi,vj,n,n1,n2)
      dimension
     *         vi(n1,n2,n),vj(n1,n2,n)
      data pi/3.1415926/
      np=n1-1
      dfi=pi*2.0/n2
      dtet=pi/np
      n6=n2/2
      cosin=cos(dtet)
      do 1 k=1,n
       do 2 j=1,n6
        j1=j+n6
        fi=(j-1)*dfi
        fi180=fi+pi
c    . . . north pole
        as=sin(cosin*fi)
        ac=cos(cosin*fi)
        asin pi=sin(cosin*fi180)
        acos pi=cos(cosin*fi180)
        a1=vj(2,j,k)*ac+
     *     vi(2,j,k)*as
        a2=vj(2,j1,k)*acos pi+
     *     vi(2,j1,k)*asin pi
        a=0.5*(a1+a2)
        b1=vj(2,j,k)*as-vi(2,j,k)*ac
        b2=vj(2,j1,k)*asin pi-
     *     vi(2,j1,k)*acos pi
        b=(b1+b2)*0.5
        vi(1,j,k)=a*sin(fi)-b*cos(fi)
        vi(1,j1,k)=-vi(1,j,k)
        vj(1,j,k)=a*cos(fi)+b*sin(fi)
        vj(1,j1,k)=-vj(1,j,k)
c    . . . south pole
        a1=vj(np,j,k)*ac-
     *     vi(np,j,k)*as
        a2=vj(np,j1,k)*acos pi-
     *     vi(np,j1,k)*asin pi
        a=0.5*(a1+a2)
        b1=-vj(np,j,k)*as-vi(np,j,k)*ac
        b2=-vj(np,j1,k)*asin pi-
     *     vi(np,j1,k)*acos pi
        b=(b1+b2)*0.5
        vi(n1,j,k)=-a*sin(fi)-b*cos(fi)
        vi(n1,j1,k)=-vi(n1,j,k)
        vj(n1,j,k)=a*cos(fi)-b*sin(fi)
        vj(n1,j1,k)=-vj(n1,j,k)
   2   continue
   1  continue
      return
      end
