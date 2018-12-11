c
c    text  <  ---  >   integer
c
      subroutine kodir(text,ntext,k,texto,ntexto)
      integer n(80),ntext(80),texto(80)
      character*1 text(80),ntexto(80)
      character*1 nn(80)/'a','b','c','d','e','f','g','h',
     *                 'i','j','k','l','m','n','o','p',
     *                 'q','r','s','t','u','v','w','x',
     *                 'y','z','1','2','3','4','5','6',
     *                 '7','8','9','0','+','|','"','#',
     *                 '.','%',',','''','(',')',' ','=',
     *                 ';','-',':','*','>','<','Z','Y',
     *                 'X','W','V','U','T','S','R','Q',
     *                 'P','O','N','M','L','K','J','I',
     *                 'H','G','F','E','D','C','B','A'/
c     character*1 nn(80)/'a','b','c','d','e','f','g','h',
c    *                 'i','j','k','l','m','n','o','p',
c    *                 'q','r','s','t','u','v','w','x',
c    *                 'y','z','1','2','3','4','5','6',
c    *                 '7','8','9','0','+','|','U','#',
c    *                 'X','Y','Z','''','(',')',' ','=',
c    *                 ';','-',':','*','>','<','?',',',
c    *                 '.','/','W','V','A','B','C','D',
c    *                 'E','F','G','H','I','J','K','L',
c    *                 'M','N','O','P','Q','R','S','T'/
      if(k.eq.1) go to 10
c*** text -> integer    ;    k=0
      do 3 i=1,80
        do 4 j=1,80
           if(text(i).ne.nn(j))go to 5
              ntext(i)=j
              go to 3
    5      continue
    4   continue
    3 continue
      go to 9
c*** integer -> text k=1
   10 continue
      do 13 i=1,80
        do 14 j=1,80
           if(texto(i).ne.j)go to 15
              ntexto(i)=nn(j)
              go to 17
   15      continue
   14   continue
   17 continue
       if(texto(i).lt.1.or.texto(i).gt.80)ntexto(i)=nn(47)
   13 continue
    9 continue
      return
      end

