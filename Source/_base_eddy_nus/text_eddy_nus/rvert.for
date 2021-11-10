      subroutine rvert(nf,kpart,ddolgt,ntsl,nl,kdf,ldor,isp,
     *       pole,nv,idt,vert,mast)
      dimension pole(ldor/4),kdf(20),vert(kpart,nv)
     *       ,ntsl(nl),mast(40)
      allocatable pp(:,:)

	logical readfl
      
	allocate (pp(kpart,nv))

	npp=kpart*nv
      readfl=.true.
      md=0
      nomsl=1
      dolm=0.
 !     do 8 i=1,nv
 !       do 9 j=1,kpart
          vert=0.
 !   9   continue
 !   8 continue
      do 10 k=1,idt
        call wwt(readfl,nf,kpart,dolm,ddolgt,nomsl,ntsl,nl,
     *           kdf,ldor,isp,md,pp,pole,npp,mast)
        do 11 i=1,nv
          do 12 j=1,kpart
            vert(j,i)=vert(j,i)+pp(j,i)
   12     continue
   11   continue
        dolm=dolm+ddolgt
   10 continue
      do 13 i=1,nv
        do 14 j=1,kpart
          vert(j,i)=vert(j,i)/float(idt)
   14   continue
   13 continue
      deallocate (pp)
      return
      end
