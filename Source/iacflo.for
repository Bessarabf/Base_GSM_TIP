c     . . . ��⮪� � ��. �ࣨ� ���஭��� ���믠���
      subroutine iacflo(utn,sole,solu,bmpz,bmpy,mas,
     *                  csol,vsol,fa0,pkp0,ap0,ae0,dst0,al0,au0,
     *                  fa,pkp,ap,ae,dst,al,au,solen,nsu,nse,ps)
      dimension sole(nse),solu(nsu),mas(10),solen(nse),ps(10)
c . . .��� ���஭�  70 �ࠤ
      ps(1 )=1.e8
c      ps(1 )=5.e8
c     ps(1 )=5.e9
c     ps(1 )=2.5e9
c     ps(1 )=1.5e9
      ps(2 )=1.
      ps(3 )=0.20e3
c . . .���⪨� ���஭� 70 �ࠤ
      ps(4 )=1.e8
c     ps(4 )=5.e8
c     ps(4 )=4.e8
      ps(5 )=1.
      ps(6 )=3.e3
c     ps(6 )=5.e3
c . . .��� ���஭�  80 �ࠤ
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
c . . . ������� ������������ ���������� ...
      if(mas(2).eq.0) then
c       . . .������� �� �室��� ���祭��
        fa=fa0
        pkp=pkp0
        ap=ap0
        ae=ae0
        dst=dst0
        al=al0
        au=au0
c       . . .�����뢠����
      else
        fa=1
        pkp=1
        ap=1
        ae=1
        dst=1
        al=1
        au=1
      end if
c    . . . �������� ��������� � ����������� ����
c    . . . vtau0 - �६� ��砫� ��������� ����⠭�� (cek)
c    . . . vtau1 - ��砫� �ᯮ���樠�쭮�� �������
c    . . . alfp  - ���୮� �६� 㬥��襭��
c    . . . pmax  - ���ᨬ��쭮� 㢥��祭�� ��ࢮ��砫쭮�� ��⮪� (ࠧ)
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
