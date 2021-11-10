! Расчитывается полная энергия, передаваемая сверхтепловыми электронами
!  тепловой плазме на плазмосферном участке трубки
!
	subroutine ctoen(aqfn,aqfs,ns,pga,pun,wm,qp)
	integer ns
	real aqfn(*),aqfs(*),pga,pun,wm(*),qp
	integer i,k
	real a,b,f,p,w,ww,x(23),y(23)
	do k=1,ns
	  i=ns-k+1
	  w=wm(k)
	  ww=w*w
	  p=ww-pun
	  if(p.lt.0.) p=0.
	  a=pga/ww
	  b=w-sqrt(p)*cpro(a)
	  f=aqfn(k)+aqfs(k)
	  y(i)=b*f
	  x(i)=w
	end do
	qp=unta(1,ns,x,y)
	end