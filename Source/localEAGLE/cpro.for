	function cpro(g)
	if(g.lt.1.e-5) then
	  a=1.
	elseif(g.gt.1.e15) then
	  a=0.
	else
        a=((exp(-g)-1.)/g+1.)*2./g
	end if
	cpro=a
	end
