c     . . . потоки и хар. энергия электронных высыпаний
      subroutine iacflo(utn,sole,solu,bmpz,bmpy,mas,
     *                  csol,vsol,fa0,pkp0,ap0,ae0,dst0,al0,au0,
     *                  fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps)
      dimension sole(nse),solu(nsu),mas(10),solen(nse),ps(10)
c . . .мягкие электроны  70 град
      ps(1 )=1.e8
c      ps(1 )=5.e8
c     ps(1 )=5.e9
c     ps(1 )=2.5e9
c     ps(1 )=1.5e9
      ps(2 )=1.
      ps(3 )=0.20e3
c . . .жесткие электроны 70 град
      ps(4 )=1.e8
c     ps(4 )=5.e8
c     ps(4 )=4.e8
      ps(5 )=1.
      ps(6 )=3.e3
c     ps(6 )=5.e3
c . . .мягкие электроны  80 град
	ps(7 )=1.e8
C      ps(7 )=5.e8
c     ps(7 )=1.e9
c     ps(7 )=2.e9
c     ps(7 )=4.e9
      ps(8 )=0.1e3
c     ps(8 )=0.05e3
      ps(9 )=0.e8
c     ps(9 )=3.e8
      ps(10)=0.05e3
c . . . ИНДЕКСЫ ГЕОМАГНИТНОЙ АКТИВНОСТИ ...
      if(mas(2).eq.0) then
c       . . .берутся из входных значений
        fa=fa0
        pkp=pkp0
        ap=ap0
        ae=ae0
        dst=dst0
        al=al0
        au=au0
c       . . .рассчитываются
      else
        fa=1
        pkp=1
        ap=1
        ae=1
        dst=1
        al=1
        au=1
      end if
c    . . . УСИЛЕНИЕ ВЫСЫПАНИЙ И ПОСТЕПЕННЫЙ СПАД
c    . . . vtau0 - время начала линейного нарастания (cek)
c    . . . vtau1 - начало экспоненциального падения
c    . . . alfp  - харктерное время уменьшения
c    . . . pmax  - максимальное увеличение первоначального потока (раз)
c     vtau0=57000.
c     if(utn.ge.vtau0) then
c        vtau1=vtau0+600.
c        alfp=2400.
c        pmax=29.
c        ps(4)=ps(4)*corp(vtau0,vtau1,utn,alfp,pmax)
c        ps(6)=5.e3
c        write(*,*) 'corp'
c     endif
      return
      end
