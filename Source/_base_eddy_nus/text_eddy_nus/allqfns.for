!	���பᨬ��� �� ���䮭��� �ࣥ��᪮�� ᯥ��� ��⮪�� 
!   ᢥ��⥯����� ���஭�� � ����� ��������, ���ࠢ������ �� ��������
!   � ����������.
!
      subroutine allqfns(ns,qfn,qfs,wm,aqfn,aqfs)
      dimension aqfn(*),aqfs(*),wm(*)
      data eo/10./,ea/10./
      an=qfn/ea
      as=qfs/ea
      do kw=1,ns
        w=wm(kw)
        b=w/eo
        exwe=b*exp(-b)
        aqfn(kw)=an*exwe
        aqfs(kw)=as*exwe
      end do
      return
      end