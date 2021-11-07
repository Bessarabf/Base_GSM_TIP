      DOUBLE PRECISION FUNCTION DERF(Y)
C
      DOUBLE PRECISION Y,X,RES,XSQ,XNUM,XDEN,XI
      DIMENSION  P(3),Q(2),P1(5),Q1(4),P2(3),Q2(2)
      DATA P(1)/.316652891D0/,P(2)/1.72227577D0/,P(3)/21.385332D0/
      DATA         Q(1)/7.84374571D0/,Q(2)/18.9522572D0/
      DATA P1(1)/.563169619D0/,P1(2)/3.03179934D0/,P1(3)/6.86501848D0/,
     *             P1(4)/7.37388831D0/,P1(5)/4.3187787D-5/
      DATA         Q1(1)/5.35421679D0/,Q1(2)/12.7955295D0/,
     *             Q1(3)/15.1849082D0/,Q1(4)/7.3739609D0/
      DATA         P2(1)/-5.16882262D-2/,P2(2)/-.196068974D0/,
     *             P2(3)/-4.25799644D-2/
      DATA         Q2(1)/.921452412D0/,Q2(2)/.150942071D0/
      DATA         XMIN/1.0D-5/,XLARGE/5.6D0/,XBIG/9.25D0/,
     *             SQRPI/.564189584D0/
      KRET=0
    5 X=Y
   10 ISW=1
      IF(X.GE.0.0D0) GO TO 15
      ISW=-1
      X=-X
   15 IF(X.LT. .477D0) GO TO 20
      IF(X.LE.4.0D0) GO TO 35
      IF(KRET*ISW.GT.0) GO TO 45
      IF(X.LT.XLARGE) GO TO 50
      RES=1.0D0
      IF(KRET.EQ.0) GO TO 65
      RES=2.0D0
      GO TO 70
C          ABS(Y) .LT. .477
   20 IF(X.LT.XMIN) GO TO 25
      XSQ=X*X
      XNUM=(P(1)*XSQ+P(2))*XSQ+P(3)
      XDEN=(XSQ+Q(1))*XSQ+Q(2)
      RES=X*XNUM/XDEN
      GO TO 30
   25 RES=X*P(3)/Q(2)
   30 IF(ISW.EQ.-1) RES=-RES
      IF(KRET.EQ.0) GO TO 80
      RES=1.0D0 - RES
      GO TO 70
C          .477 .LE. ABS(Y) .LE. 4.0D0
   35 XSQ=X*X
      XNUM=P1(5)*X+P1(1)
      XDEN=X+Q1(1)
      DO 40 I=2,4
         XNUM=XNUM*X+P1(I)
         XDEN=XDEN*X+Q1(I)
   40 CONTINUE
      RES=XNUM/XDEN
      GO TO 55
C            4.0 .LT. ABS(Y)
   45 IF(X.GT.XBIG) GO TO 75
   50 XSQ=X*X
      XI=1.0D0/XSQ
      XNUM=((P2(1)*XI+P2(2))*XI+P2(3))
      XDEN=(XI+Q2(1))*XI+Q2(2)
      RES=(SQRPI+XI*XNUM/XDEN)/X
   55 RES=RES*DEXP(-XSQ)
      IF(KRET.EQ.0) GO TO 60
      IF(ISW.EQ.-1) RES=2.0D0 - RES
      GO TO 70
   60 RES=1.0D0 - RES
   65 IF(ISW.EQ.-1) RES=-RES
      GO TO 80
   70 IF(KRET.EQ.2) RES=RES*0.5D0
      GO TO 80
   75 RES=0.0D0
   80 DERF=RES
      RETURN
      END

