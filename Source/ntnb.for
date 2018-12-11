      subroutine ntnb(nx,ht,tt,cio,cih,cihe,cim,b1,snt,snb)

      dimension ht(*),tt(*),cio(*),cih(*),cihe(*),cim(*)
      is=18
      in=nx-18
      snb=0.
      snt=0.
      hm=ht(is)
      tm=tt(is)
      bm=bdip(hm,tm)
      b1=bm
      cem=cio(is)+cih(is)+cihe(is)
      if(hm.le.1.e8)cem=cem+cim(is)
      do 1 i=is,in
        ip=i+1
        hp=ht(ip)
        alt=(hm+hp)*.5
        tp=tt(ip)
        dst=ds(hm,alt,hp,tm,tp)*.5
        bp=bdip(hp,tp)
        cep=cio(ip)+cih(ip)+cihe(ip)
        if(hp.le.1.e8)cep=cep+cim(ip)
        snt=snt+(cem+cep)*dst
        snb=snb+(cem/bm+cep/bp)*dst
        hm=hp
        tm=tp
        bm=bp
        cem=cep
    1 continue
      snb=snb*b1
      return
      end

