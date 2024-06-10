      PROGRAM OMEGA2
C------THIS PROGRAM FINDS THE TWIN LAWS IN A MANNER SIMILAR TO
C------TO THAT DESCRIBED IN YVON LE PAGE'S ORIGINAL 1982 PAPER
C------ON CELL REDUCTIONS.  tHE REULT IS SIMILAR TO THAT PROPOSED
C------FOR THE OBLIQUE PROGRAM OF LE PAGE, BUT (WE BELIEVE) MORE
C------COMPLETE FOR THE SPECIFIED SCALAR PRODUCTS (VERBOSE OUTPUT)
C------OR TWIN INDEX m, (PRODUCED IN VERBOSE OR NON-VERBOSE MODE).
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOM/ AA(3,3),AINV(3,3),NH(3,1000000),HH(3,1000000),
     1AANG(1000000),TRANS(3,3),IMAX,IMAX1,IMAX2,N2
      COMMON /TFORM/ TIND(3,3),RIND(3,3)
      COMMON /NEAT/ JD(3,500,64),JR(3,500,64),JDT(3,500,64),JRT(3,500,
     164),JMULT(500,64),JTWIN(500,64),OBL(500,64),JC(64)
      DIMENSION CELL(3), CELANG(3), IPT(3), IHXT(3,1000000)
      DIMENSION P(3), IP(3), HX(3,1000000), IHX(3,1000000), IH(3), HI(3)
     1, AIPT(3), AIHXT(3), DUM(3,3), ARG(3,3), BRG(3,3), CRG(3,3)
      COMMON /FINAL/ JDN(3,500,32),JRN(3,500,32),JDTN(3,500,32),JRTN(3,
     1500,32),JMULTN(500,32),JTWINN(500,32),OBLN(500,32),JCNEW(32),
     2NTEST(500),NLAWS
      CHARACTER ANS*1,DATE*8,TIME*8,SITE*6,ATMOD*1,DUMA*1,TEST*4,TITLE*
     132,FILOUT*32
      DATA IOUT,IOUT2,INP,NVERBOSE,Z1/8,15,5,0,1.0D00/
      INP=4
      IIN=4
      WRITE (*,10)
   10 FORMAT (//' Program OMEGA, adapted from and inspired by ','CREDUC8
     11: ',/,' Yvon Le Page, J. Appl. Cryst. (1982) 15, ','255-259,'/' V
     2ersion 2.50 by Bruce Foxman, Brandeis University, ','25 MAY 2020')
   20 FORMAT (A1)
      WRITE (*,30) 
   30 FORMAT (//' Enter file name for output [DEFAULT=OMEGA.OUT]: ',$)
      READ (*,'(A32)') FILOUT
      IF (FILOUT.EQ.'') FILOUT='OMEGA.OUT'
      OPEN (8,STATUS='UNKNOWN',FILE=FILOUT,ACCESS='SEQUENTIAL',FORM='FOR
     1MATTED')
c       WRITE(*,58)
C 58     FORMAT(/' Omit verbose outout in file? [y or n; default = n]: '
C     1 ,$)
C       READ(*,'(A4)') TEST
      IF (TEST.EQ.''.OR.TEST.EQ.'N'.OR.TEST.EQ.'n') GO TO 40
      NVERBOSE=1
   40 ANGMAX=3.0D00
      IRET=1
      NEWOLD=1
      NLAWS=0
      DO 50 KK=1,64
   50   JC(KK)=0
      WRITE (6,60)
   60 FORMAT (/' Do you want to change the maximum acceptable obliquity'
     1,' (omega) ?'/,' (Default value is 3.0) (Y/N): ',$)
      READ (5,20,end=120) ans
      IF ((ans.ne.'N').and.(ans.ne.'n').and.(ans.ne.'Y').and.(ans.ne.'y'
     1)) GO TO 40
      IF ((ans.eq.'Y').OR.(ans.eq.'y')) THEN
        WRITE (6,70)
   70   FORMAT (/' What is the largest maximum omega(degrees) to conside
     1r?',' ',$)
        READ (5,*,end=120) angmax
        IF (angmax.le.0.0D00) angmax=3.0D00
      END IF
      WRITE (6,80)
   80 FORMAT (/' Do you want to search for twin indices (m) higher than'
     1,' 2 ?',/,' (y/n, default = n) ? ',$)
      READ (5,20,end=120) ans
      IF ((ans.eq.'Y').or.(ans.eq.'y')) THEN
        WRITE (6,90)
   90   FORMAT (/' What is the largest value of DOT (the scalar product)
     1',' to use? ',//,' Note that the maximum reported twin index, m, '
     2   /' will be one-','half the value of DOT (or DOT-1 if DOT is odd
     3):  ',$)
        READ (5,*,end=120) imax
        WRITE (6,100)
  100   FORMAT (//' See http://reference.iucr.org/dictionary/','Twin_ind
     1ex',/' for further information.'/)
        IF (imax.le.1) imax=2
      ELSE
        IMAX=2
      END IF
      IMAX1=IMAX+1
      IMAX2=2*IMAX+1
      IRET=0
      WRITE (IOUT,10) 
      WRITE (iout,110) angmax,imax
      WRITE (*,110) angmax,imax
  110 FORMAT (/' Maximum obliquity = ',f6.2,/,' Maximum scalar ','produc
     1t (DOT) = ',i3)
      CALL CINPUT (INP, IOUT, CELL, CELANG, TRANS, IRET, ATMOD, ITEST)
  120 CONTINUE
      IF (IRET.EQ.1) THEN
        WRITE (6,130)
  130   FORMAT (//' La reduction d''Yvon Le Page est terminee')
        STOP
      END IF
      IF (NVERBOSE.EQ.0) WRITE (IOUT,140)
  140 FORMAT (//11X,'TRADITIONAL CREDUC81-STYLE OUTPUT FOLLOWS:'/)
      IF (ITEST.EQ.0) GO TO 170
      IF (NVERBOSE.EQ.0) WRITE (IOUT,150)
  150 FORMAT (/15X,'ROWS',22X,'PRODUCTS'//11X,'BUERGER CELL',37X,'     I
     1NPUT CELL')
      IF (NVERBOSE.EQ.0) WRITE (IOUT,160)
  160 FORMAT (/,7X,'DIRECT',7X,'RECIPROCAL',9X,'DOT',4X,' OMEGA',8X,'DIR
     1ECT',8X,'RECIPROCAL')
      GO TO 200
  170 CONTINUE
      IF (NVERBOSE.EQ.0) WRITE (IOUT,180)
  180 FORMAT (/15X,'ROWS',22X,'PRODUCTS'//11X,'BUERGER CELL')
      IF (NVERBOSE.EQ.0) WRITE (IOUT,190)
  190 FORMAT (/,7X,'DIRECT',7X,'RECIPROCAL',9X,'DOT',4X,' OMEGA')
C------GET THE ORTHOGONAL REPRESENTATION OF THE PRIMITIVE UNIT VECTORS
  200 DO 210 J=1,3
        ARG(J,1)=CELL(J)
        BRG(J,1)=CELANG(J)
  210 CONTINUE
      CALL MATRIX (ARG, BRG, AA, DUM, 'ORTHOG')
      CALL MATRIX (AA, AINV, DUM, DUM, 'INVERT')
C------GENERATE ROW MILLER INDICES NO LARGER THAN IMAX
      ITOT=0
      IFL=0
      N2=0
C------GENERATE ALL INDEX COMBINATIONS
      DO 230 IA=1,IMAX1
        IH(1)=IA-1
        DO 230 IB=1,IMAX2
        IH(2)=IB-IMAX-1
        DO 230 IC=1,IMAX2
        IH(3)=IC-IMAX-1
C------REJECT THE ONES THAT ARE NOT COPRIME
        CALL COPRIM (IH, KA, IMAX)
        IF (KA.NE.0) GO TO 230
C        IF((IABS(IH(1))+IABS(IH(2))+IABS(IH(3))).EQ.5.AND.IMAX.EQ.2)
C     1  GO TO 20
        ITOT=ITOT+1
C------GET THE ORTHOGONAL COMPONENTS OF THE ROW
        HI(1)=IH(1)
        HI(2)=IH(2)
        HI(3)=IH(3)
        BRG(1,1)=HI(1)
        BRG(2,1)=HI(2)
        BRG(3,1)=HI(3)
        CALL MATRIX (AA, BRG, HX(1,ITOT), DUM, 'MVMULT')
        DO 220 JA=1,3
  220     IHX(JA,ITOT)=IH(JA)
        IF (ITOT.LT.1000000) GO TO 230
        GO TO 240
  230 CONTINUE
C------IF WE PASS HERE WE HAVE EXHAUSTED THE POSSIBLE ROWS
      IFL=1
C------GENERATE TRICLINIC MILLER INDICES FOR PLANES
C        ITEST=1
  240 DO 360 I=-imax1,IMAX1
        IP(1)=I-1
        DO 360 J=1,IMAX2
        IP(2)=J-IMAX-1
        DO 360 K=1,IMAX2
        IP(3)=K-IMAX-1
        CALL COPRIM (IP, KA, IMAX)
C        IF((IABS(IP(1))+IABS(IP(2))+IABS(IP(3))).EQ.5.AND.IMAX.EQ.2)
C     1  GO TO 10
        IF (KA.NE.0) GO TO 360
C------GET THE ORTHOGONAL COMPONENTS OF THE NORMAL TO THE PLANE
        HI(1)=IP(1)
        HI(2)=IP(2)
        HI(3)=IP(3)
        ARG(1,1)=HI(1)
        ARG(2,1)=HI(2)
        ARG(3,1)=HI(3)
        CALL MATRIX (ARG, AINV, CRG, DUM, 'VMMULT')
        P(1)=CRG(1,1)
        P(2)=CRG(2,1)
        P(3)=CRG(3,1)
        IF (ITOT.EQ.0) GO TO 430
C------GET ONE ROW
        DO 350 L=1,ITOT
C------CALCULATE THE MULTIPLICITY OF THE CELL DEFINED BY THE MESH ON THE
C------PLANE AND THE TRANSLATION ALONG THE ROW
          MULT=IABS(IHX(1,L)*IP(1)+IHX(2,L)*IP(2)+IHX(3,L)*IP(3))
          IF (MULT.EQ.0.OR.MULT.GT.IMAX) GO TO 350
C------CALCULATE THE ANGLE BETWEEN THE ROW AND THE NORMAL TO THE PLANE
          CALL DOTP (P(1), P(2), P(3), HX(1,L), HX(2,L), HX(3,L), ANG)
C        IF(ANG.GT.ANG1TEST.AND.ANG.LT.180.) ITEST=-1
C        IF(ITEST.NE.-1) GO TO 412
C        GO TO 413
          IF (ANG.GT.ANGMAX) GO TO 350
c        ANG=ATAN(SQRT(ANG))*180/3.1415927
C------THIS IS A POTENTIAL EVEN-ORDER AXIS, WE PRINT IT AND SAVE IT
C        WRITE(IOUT,11)(IHX(JA,L),JA=1,3),IP,MULT,ANG
  250     FORMAT (2X,'[',3I4,']',2X,'(',3I4,')',I10,F10.4)
          IF (ITEST.EQ.0) GO TO 310
          IFAC=1
          IPT(1)=RIND(1,1)*IP(1)+RIND(1,2)*IP(2)+RIND(1,3)*IP(3)
          IPT(2)=RIND(2,1)*IP(1)+RIND(2,2)*IP(2)+RIND(2,3)*IP(3)
          IPT(3)=RIND(3,1)*IP(1)+RIND(3,2)*IP(2)+RIND(3,3)*IP(3)
          AIHXT(1)=(TIND(1,1)*IHX(1,L)+TIND(2,1)*IHX(2,L)+TIND(3,1)*
     1     IHX(3,L))
          AIHXT(2)=(TIND(1,2)*IHX(1,L)+TIND(2,2)*IHX(2,L)+TIND(3,2)*
     1     IHX(3,L))
          AIHXT(3)=(TIND(1,3)*IHX(1,L)+TIND(2,3)*IHX(2,L)+TIND(3,3)*
     1     IHX(3,L))
          DO 260 JK=1,3
  260       AIPT(JK)=IPT(JK)
          IF (DMOD(AIHXT(1),Z1).NE.O.OR.DMOD(AIHXT(2),Z1)
     1     .NE.0.OR.DMOD(AIHXT(3),Z1).NE.0) IFAC=2
          DO 270 JJ=1,3
  270       AIHXT(JJ)=IFAC*AIHXT(JJ)
          CALL COPRIME2 (AIHXT, AIPT)
          DO 280 JK=1,3
            IHXT(JK,L)=AIHXT(JK)
  280       IPT(JK)=AIPT(JK)
          WRITE (IOUT,300) (IHX(JA,L),JA=1,3),IP,MULT,ANG,(IHXT(JA,L),
     1     JA=1,3),IPT
          JC(MULT)=JC(MULT)+1
          DO 290 JJ=1,3
            JDT(JJ,JC(MULT),MULT)=IHXT(JJ,L)
            JRT(JJ,JC(MULT),MULT)=IPT(JJ)
            JD(JJ,JC(MULT),MULT)=IHX(JJ,L)
  290       JR(JJ,JC(MULT),MULT)=IP(JJ)
          JMULT(JC(MULT),MULT)=MULT
          JTWIN(JC(MULT),MULT)=MULT
          IF (MOD(MULT,2).EQ.0) JTWIN(JC(MULT),MULT)=JTWIN(JC(MULT),
     1     MULT)/2
          OBL(JC(MULT),MULT)=ANG
          GO TO 330
  300     FORMAT (2X,'[',3I4,']',2X,'(',3I4,')',I10,F10.4,4X,'[',3I4,']'
     1     ,2x,'(',3I4,')')
  310     IF (NVERBOSE.EQ.0) WRITE (IOUT,250) (IHX(JA,L),JA=1,3),IP,
     1     MULT,ANG
          JC(MULT)=JC(MULT)+1
          DO 320 JJ=1,3
            JD(JJ,JC(MULT),MULT)=IHX(JJ,L)
  320       JR(JJ,JC(MULT),MULT)=IP(JJ)
          JMULT(JC(MULT),MULT)=MULT
          JTWIN(JC(MULT),MULT)=MULT
          IF (MOD(MULT,2).EQ.0) JTWIN(JC(MULT),MULT)=JTWIN(JC(MULT),
     1     MULT)/2
          OBL(JC(MULT),MULT)=ANG
  330     HMOD=SQRT(HX(1,L)**2+HX(2,L)**2+HX(3,L)**2)
          N2=N2+1
          DO 340 NX=1,3
            NH(NX,N2)=IHX(NX,L)
            HH(NX,N2)=HX(NX,L)/HMOD
  340       AANG(N2)=ANG
  350   CONTINUE
  360 CONTINUE
      IF (NVERBOSE.EQ.0) WRITE (IOUT,370)
  370 FORMAT (//11X,'TRADITIONAL CREDUC81-STYLE OUTPUT, SORTED ON DOT,',
     1' FOLLOWS: '/)
      IF (ITEST.EQ.0) GO TO 380
      IF (NVERBOSE.EQ.0) WRITE (IOUT,150)
      IF (NVERBOSE.EQ.0) WRITE (IOUT,160)
      GO TO 390
  380 IF (NVERBOSE.EQ.0) WRITE (IOUT,180)
      IF (NVERBOSE.EQ.0) WRITE (IOUT,190)
  390 DO 420 J=1,IMAX
        CALL SORT (JC(J), J, ATMOD, NEWOLD, ITEST)
        DO 410 K=1,JC(J)
          IF (ITEST.EQ.0) GO TO 400
          IF (NVERBOSE.EQ.0) WRITE (IOUT,300) (JD(L,K,J),L=1,3),(JR(L,K,
     1     J),L=1,3),JMULT(K,J),OBL(K,J),(JDT(L,K,J),L=1,3),(JRT(L,K,J),
     2     L=1,3)
          GO TO 410
  400     IF (NVERBOSE.EQ.0) WRITE (IOUT,250) (JD(L,K,J),L=1,3),(JR(L,K,
     1     J),L=1,3),JMULT(K,J),OBL(K,J)
  
  410   CONTINUE
  420 CONTINUE
      CALL REGROUP (IMAX, ATMOD, IOUT, ITEST)
      IF (IFL.EQ.1) GO TO 430
C------THERE ARE MORE ROWS TO BE GENERATED
      ITOT=0
      GO TO 230
C------ORDER THE ROWS ON THE ANGLE WITH THE NORMAL TO THE PLANE
  430 IF (N2.LT.2) GO TO 470
      DO 460 I=1,N2-1
        ANMAX=AANG(I)
        MAX=I
        DO 440 J=I+1,N2
          IF (AANG(J).GE.ANMAX) GO TO 440
          ANMAX=AANG(J)
          MAX=J
  440   CONTINUE
        DO 450 J=1,3
          ISAVE=NH(J,I)
          NH(J,I)=NH(J,MAX)
          NH(J,MAX)=ISAVE
          SAVE=HH(J,I)
          HH(J,I)=HH(J,MAX)
  450     HH(J,MAX)=SAVE
        AANG(MAX)=AANG(I)
        AANG(I)=ANMAX
  460 CONTINUE
C------FIND THE CRYSTAL SYSTEM
  470 WRITE (IOUT,480)
  480 FORMAT (/)
C        CALL FNDSYS(IOUT,HH,N2,N2)
C        WRITE(IOUT,566)
C566     FORMAT(//)
      WRITE (6,490)
  490 FORMAT (/' New Run? (Y/N) : ',$)
      READ (5,*,end=120) ans
      IF ((ans.eq.'n').or.(ans.eq.'N')) WRITE (*,500)
  500 FORMAT (//' END OF SESSION'/)
      IF ((ans.eq.'n').or.(ans.eq.'N')) STOP
      GO TO 40
      END
      SUBROUTINE MATRIX (A, B, C, D, WHAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 WHAT
      DIMENSION A(3,3), B(3,3), C(3,3), D(3,3), E(3,3), V(3)
      IF (WHAT.EQ.'INVERT') GO TO 20
      IF (WHAT.EQ.'MATMUL') GO TO 40
      IF (WHAT.EQ.'MATVEC') GO TO 70
      IF (WHAT.EQ.'VECMAT') GO TO 100
      IF (WHAT.EQ.'SCALPR') GO TO 130
      IF (WHAT.EQ.'LENGTH') GO TO 150
      IF (WHAT.EQ.'ORTHOG') GO TO 170
      IF (WHAT.EQ.'DETERM') GO TO 190
      IF (WHAT.EQ.'MVMULT') GO TO 210
      IF (WHAT.EQ.'VMMULT') GO TO 230
      IF (WHAT.EQ.'TRNSPS') GO TO 250
      WRITE (4,10) WHAT
   10 FORMAT (' MATRIX OPERATION ',A6,' IS NOT PROGRAMMED')
C------INVERT 3X3 MATRIX A, PUT THE RESULT IN B
   20 E(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      E(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      E(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      E(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      E(2,2)=A(1,1)*A(3,3)-A(1,3)*A(3,1)
      E(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      E(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
      E(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      E(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      DMAT=A(1,1)*E(1,1)+A(1,2)*E(2,1)+A(1,3)*E(3,1)
      DO 30 I=1,3
        DO 30 J=1,3
   30   B(I,J)=E(I,J)/DMAT
      GO TO 280
C------MULTIPLY 3X3 MATRICES A AND B, STORE RESULT IN C
   40 DO 50 I=1,3
        DO 50 J=1,3
        E(I,J)=0.0
        DO 50 K=1,3
   50   E(I,J)=E(I,J)+A(I,K)*B(K,J)
      DO 60 I=1,3
        DO 60 J=1,3
   60   C(I,J)=E(I,J)
      GO TO 280
C------MULTIPLY MATRIX A BY VECTOR B, STORE DIRECTION COSINES OF RESULT
   70 DO 80 I=1,3
        V(I)=0.
        DO 80 J=1,3
   80   V(I)=V(I)+A(I,J)*B(J,1)
      VMOD=SQRT(V(1)**2+V(2)**2+V(3)**2)
      DO 90 I=1,3
   90   C(I,1)=V(I)/VMOD
      GO TO 280
C------MULTIPLY VECTOR A BY MATRIX B, STORE DIRECTION COSINES OF RESULT
  100 DO 110 I=1,3
        V(I)=0.
        DO 110 J=1,3
  110   V(I)=V(I)+B(J,I)*A(J,1)
      VMOD=SQRT(V(1)**2+V(2)**2+V(3)**2)
      DO 120 I=1,3
  120   C(I,1)=V(I)/VMOD
      GO TO 280
C------SCALAR PRODUCT OF VECTORS A AND B
  130 S=0
      DO 140 I=1,3
  140   S=S+A(I,1)*B(I,1)
      C(1,1)=S
      GO TO 280
C------LENGTH OF VECTOR B WHEN A IS THE METRIC MATRIX
  150 DO 160 I=1,3
        V(I)=0.
        DO 160 J=1,3
  160   V(I)=V(I)+A(I,J)*B(J,1)
      C(1,1)=SQRT(V(1)**2+V(2)**2+V(3)**2)
      GO TO 280
C------GET THE METRIC MATRIX C CORRESPONDING TO CELL EDGES A AND ANGLES
  170 COSGAS=(DCOSD(B(1,1))*DCOSD(B(2,1))-DCOSD(B(3,1)))
      COSGAS=COSGAS/(DSIND(B(1,1))*DSIND(B(2,1)))
      SINGAS=SQRT(1.-COSGAS**2)
      E(1,1)=A(1,1)*DSIND(B(2,1))*SINGAS
      E(1,2)=0
      E(1,3)=0
      E(2,1)=-A(1,1)*DSIND(B(2,1))*COSGAS
      E(2,2)=A(2,1)*DSIND(B(1,1))
      E(2,3)=0
      E(3,1)=A(1,1)*DCOSD(B(2,1))
      E(3,2)=A(2,1)*DCOSD(B(1,1))
      E(3,3)=A(3,1)
      DO 180 I=1,3
        DO 180 J=1,3
  180   C(I,J)=E(I,J)
      GO TO 280
C------CALCULATE THE DETERMINANT D OF THE VECTORS A,B,C
  190 DET=0.
      DO 200 I=1,3
        J=I+1
        IF (J.EQ.4) J=1
        K=6-I-J
  200   DET=DET+A(I,1)*(B(J,1)*C(K,1)-B(K,1)*C(J,1))
      D(1,1)=DET
      GO TO 280
C------MULTIPLY MATRIX A BY VECTOR B, STORE RESULT IN C
  210 DO 220 I=1,3
        C(I,1)=0
        DO 220 J=1,3
  220   C(I,1)=C(I,1)+A(I,J)*B(J,1)
      GO TO 280
C------MULTIPLY VECTOR A BY MATRIX B, STORE RESULT IN C
  230 DO 240 I=1,3
        C(I,1)=0.
        DO 240 J=1,3
  240   C(I,1)=C(I,1)+A(J,1)*B(J,I)
      GO TO 280
C------TRANSPOSE A AND PUT IT IN B
  250 DO 260 I=1,3
        DO 260 J=1,3
  260   E(I,J)=A(J,I)
      DO 270 I=1,3
        DO 270 J=1,3
  270   B(I,J)=E(I,J)
      GO TO 280
  280 RETURN
      END
      SUBROUTINE BURGER (IOUT, A, ANG, IND, ITEST)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 IND(3,3)
      DIMENSION AA(3,3), A(3), ANG(3), DUM(3,3)
      COMMON /TFORM/ TIND(3,3),RIND(3,3)
      DATA INUM/0/
      ITEST=0
C   FORM THE MATRIX OF DOT PRODUCTS
      DO 10 I=1,3
        DO 10 J=1,3
        IF (I.EQ.J) AA(I,J)=A(I)*A(I)
        IF (I.NE.J) AA(I,J)=A(I)*A(J)*DCOSD(ANG(6-I-J))
   10 CONTINUE
C  LOOK FOR SHORTER TRANSLATIONS IN CELL FACES
   20 NUM=0
      DO 80 I=1,3
        DO 70 J=1,3
          IF (J.EQ.I) GO TO 70
          IS=1
          IF (AA(I,J).GT.0) IS=-1
          IS1=IS
          VMIN=0
   30     V=AA(I,J)*2*IS1+AA(J,J)*IS1**2
          IF (V.GE.VMIN) GO TO 40
          VMIN=V
          IS1=IS1+IS
          GO TO 30
C  DID WE FIND A SHORTER TRANSLATION?
   40     IS1=IS1-IS
          IF (IS1.EQ.0) GO TO 60
C  YES, WE DID. ACCEPT IT AS A CELL EDGE
          NUM=NUM+1
          INUM=INUM+1
C  TRANSFORM THE OLD-NEW INDICES
          DO 50 K=1,3
   50       IND(I,K)=IND(I,K)+IS1*IND(J,K)
C  MODIFY THE MATRIX OF DOT PRODUCTS
          AA(I,I)=AA(I,I)+AA(I,J)*2*IS1+AA(J,J)*IS1**2
          AA(I,J)=AA(I,J)+IS1*AA(J,J)
          AA(J,I)=AA(I,J)
          K=6-I-J
          AA(I,K)=AA(I,K)+IS1*AA(J,K)
          AA(K,I)=AA(I,K)
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE
C  LOOK FOR MORE TRANSFORMATIONS
      IF (NUM.GE.1) GO TO 20
C  ARE THE CROSS-TERMS OF A SAME SIGN?
   90 VAR=ABS(AA(1,2))+ABS(AA(1,3))+ABS(AA(2,3))
      IF (ABS(ABS(AA(1,2)+AA(1,3)+AA(2,3))-VAR).LE..0001*VAR) GO TO 140
C  NO, FIND THE ODD SIGN
      ISIGN=1
      IF (AA(1,2)*AA(1,3)*AA(2,3).LT.0) ISIGN=-1
C  REVERSE TWO VECTORS TO MAKE THE CELL TRIACUTE OR TRIOBTUSE
      DO 110 I=1,2
        K=I+1
        DO 100 J=K,3
          IF (AA(I,J)*ISIGN.GT.(0.)) GO TO 120
  100   CONTINUE
  110 CONTINUE
  120 K=6-I-J
C  MODIFY THE INDICES AND THE DOT PRODUCTS
      DO 130 II=1,3
        IND(I,II)=-IND(I,II)
  130   IND(J,II)=-IND(J,II)
      AA(K,J)=-AA(K,J)
      AA(J,K)=-AA(J,K)
      AA(K,I)=-AA(K,I)
      AA(I,K)=-AA(I,K)
  140 CONTINUE
C  ORDER THE DIAGONAL TERMS IN INCREASING VALUES
      INUM=0
  150 NUM=0
      DO 190 I=1,2
        IF ((AA(I,I)-AA(I+1,I+1)).LE.(0.)) GO TO 190
        NUM=NUM+1
        INUM=INUM+1
        DO 160 J=1,3
          SAVE=AA(I,J)
          AA(I,J)=AA(I+1,J)
  160     AA(I+1,J)=SAVE
        DO 170 J=1,3
          SAVE=AA(J,I)
          AA(J,I)=AA(J,I+1)
  170     AA(J,I+1)=SAVE
        DO 180 K=1,3
          SAVE=IND(I,K)
          IND(I,K)=IND(I+1,K)
  180     IND(I+1,K)=SAVE
  190 CONTINUE
      IF (NUM.NE.0) GO TO 150
C IF THE CELL IS LEFT-HANDED, REVERSE ALL AXES
      IF (MOD(INUM,2).EQ.0) GO TO 210
      DO 200 I=1,3
        DO 200 J=1,3
C IF 111 IS SHORTER THAN C, CALL IT C AND RE-REDUCE THE CELL
  200   IND(I,J)=-IND(I,J)
  210 IF (AA(1,1)+AA(2,2).GE.-2*(AA(1,2)+AA(1,3)+AA(2,3))) GO TO 230
      AA(3,3)=AA(3,3)+2*AA(3,1)+AA(1,1)
      AA(3,1)=AA(3,1)+AA(1,1)
      AA(1,3)=AA(3,1)
      AA(3,2)=AA(3,2)+AA(1,2)
      AA(2,3)=AA(3,2)
      AA(3,3)=AA(3,3)+2*AA(3,2)+AA(2,2)
      AA(3,2)=AA(3,2)+AA(2,2)
      AA(2,3)=AA(3,2)
      AA(3,1)=AA(3,1)+AA(1,2)
      AA(1,3)=AA(3,1)
      DO 220 J=1,3
  220   IND(3,J)=IND(1,J)+IND(2,J)+IND(3,J)
      GO TO 90
C  GET THE BUERGER CELL PARAMETERS
  230 RADEG=180.D00/3.14159265358979D00
      DO 240 I=1,3
  240   A(I)=SQRT(AA(I,I))
      DO 250 I=1,3
        J=MOD(I,3)+1
        K=MOD(J,3)+1
  250   ANG(I)=RADEG*DACOS(AA(J,K)/(A(J)*A(K)))
      WRITE (*,260) A,ANG
      WRITE (IOUT,260) A,ANG
  260 FORMAT (/' THE BUERGER CELL: ',6F9.4)
      WRITE (IOUT,270) ((IND(I,J),J=1,3),I=1,3)
      WRITE (*,270) ((IND(I,J),J=1,3),I=1,3)
  270 FORMAT (/' THE INPUT-TO-BUERGER CELL MATRIX'//(25X,3F6.1))
      DO 280 I=1,3
        DO 280 J=1,3
  280   TIND(I,J)=IND(I,J)
      CALL MATRIX (TIND, RIND, DUM, DUM, 'INVERT')
      WRITE (IOUT,290) ((RIND(I,J),J=1,3),I=1,3)
      WRITE (*,290) ((RIND(I,J),J=1,3),I=1,3)
  290 FORMAT (/' THE BUERGER-TO-INPUT CELL MATRIX'//(25X,3F6.1))
      DO 310 I=1,3
        DO 310 J=1,3
        IF (I.EQ.J) GO TO 300
        IF (RIND(I,J).NE.0.0) ITEST=1
        GO TO 310
  300   IF (RIND(I,J).NE.1.0) ITEST=1
  310 CONTINUE
      RETURN
      END
      SUBROUTINE CINPUT (IIN, IOUT, PRIM, ANPRIM, TRANSF, IRET, ATMOD,
     1ITEST)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(3), ALP(3), SYS(7), TRANS(3,3,7), AA(3,3), PRIM(3)
      DIMENSION ANPRIM(3), TRANSF(3,3), H(3,3), SYSL(7), DUM(3,3),
     1ARG(3,3), BRG(3,3), CRG(3,3)
      COMMON /SYST/ XSYS
      CHARACTER ATMOD*1,SYS*1,TITLE*80,SYSL*1
      DATA SYS/'P','A','B','C','I','F','R'/
      DATA SYSL/'p','a','b','c','i','f','r'/
      CHARACTER*4 XSYS
      DATA TRANS/1.,0.,0.,0.,1.,0.,0.,0.,1.,1.,0.,0.,0.,.5,.5,0.,0.,1.,
     1.5,0.,.5,0.,1.,0.,0.,0.,1.,.5,.5,0.,0.,1.,0.,0.,0.,1.,.5,.5,.5,0.,
     21.,0.,0.,0.,1.,.5,.5,0.,0.,.5,.5,.5,0.,.5,.666667,.333333,.333333,
     3-.333333,.333333,.333333,-.333333,-.666667,.333333/
      XSYS='    '
      RADEG=180.D00/3.14159265358979D00
      IF (IIN.NE.4) GO TO 150
      WRITE (6,10)
   10 FORMAT (/' Cell Parameters in Angstroms and degrees?',//,' Give in
     1 order a,b,c,alpha or cos(alpha),beta or cos(beta)',',gamma or cos
     2(gamma)',//' : ',$)
      READ (5,*,end=200) A,ALP
      CALL STANDCALC (A, ALP, IOUT)
      IF (A(1).EQ.0.) GO TO 200
c---- Also allow input of cos alpha, cos beta and cos gamma
      DO 20 i=1,3
        IF (abs(alp(i)).le.1.0) alp(i)=radeg*Dacos(alp(i))
   20 CONTINUE
   30 WRITE (6,40)
   40 FORMAT (/' Lattice mode? Give P,A,B,C,I,F or R: ',$)
      READ (5,50,end=200) ATMOD
   50 FORMAT (A1)
      DO 60 I=1,7
        IF (ATMOD.EQ.SYS(I).OR.ATMOD.EQ.SYSL(I)) GO TO 70
   60 CONTINUE
      GO TO 30
   70 CONTINUE
C      if((iin.eq.4).and.(iout.eq.6)) go to 104
      WRITE (6,80)
   80 FORMAT (//' Title for This Run? ',$)
      READ (5,190,end=200) TITLE
   90 WRITE (IOUT,100) TITLE,A,ALP,ATMOD
      WRITE (*,100) TITLE,A,ALP,ATMOD
      IF (ALP(1).EQ.90..AND.ALP(2).NE.90..AND.ALP(3).EQ.90.) XSYS='MONO'
      IF (ALP(1).NE.90..AND.ALP(2).NE.90..OR.ALP(2).NE.90..AND.ALP(3)
     1.NE.90..OR.ALP(1).NE.90..AND.ALP(3).NE.90.) XSYS='TRIC'
  100 FORMAT (//'   STRUCTURE: ',a80,//' THE INPUT CELL:',3X,6F9.4,'  (L
     1ATTICE MODE ',A1,')')
      DO 110 J=1,3
        ARG(J,1)=A(J)
        BRG(J,1)=ALP(J)
  110 CONTINUE
      CALL MATRIX (ARG, BRG, AA, DUM, 'ORTHOG')
      DO 120 N=1,3
  120   CALL MATRIX (AA, TRANS(1,N,I), H(1,N), DUM, 'MATVEC')
      DO 130 N=1,3
        CALL MATRIX (AA, TRANS(1,N,I), PRIM(N), DUM, 'LENGTH')
        J=MOD(N,3)+1
        K=6-N-J
        CALL MATRIX (H(1,J), H(1,K), CRG, DUM, 'SCALPR')
        COSNG=CRG(1,1)
  130   ANPRIM(N)=DACOS(COSNG)*RADEG
      DO 140 N=1,3
        DO 140 NN=1,3
  140   TRANSF(NN,N)=TRANS(N,NN,I)
      CALL BURGER (IOUT, PRIM, ANPRIM, TRANSF, ITEST)
      RETURN
  150 CONTINUE
      READ (IIN,*,end=200) A,ALP,ATMOD
      IF (ATMOD.EQ.' ') ATMOD='P'
      IF (A(1).EQ.0.) GO TO 200
      DO 160 i=1,7
        IF (atmod.eq.sys(i)) GO TO 170
  160 CONTINUE
      atmod=sys(1)
  170 CONTINUE
c---- Also allow input of cos alpha, cos beta and cos gamma
      DO 180 i=1,3
        IF (abs(alp(i)).le.1.0) alp(i)=radeg*Dacos(alp(i))
  180 CONTINUE
      READ (iin,190,end=200) JUNK
  190 FORMAT (A80)
      GO TO 90
  200 IRET=1
      RETURN
      END
      SUBROUTINE COPRIM (I, K, IMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION I(3)
      K=0
C        IF(I(1).NE.0)GO TO 1
C        IF(I(2).NE.0)GO TO 2
C        IF(I(3).NE.1)K=1
C2       IF(I(2).LT.0)K=1
      IF (MOD(I(1),2).NE.0) GO TO 10
      IF (MOD(I(2),2).NE.0) GO TO 10
      IF (MOD(I(3),2).NE.0) GO TO 10
      K=1
      RETURN
   10 IF (IMAX.LT.3) RETURN
      IF (MOD(I(1),3).NE.0) GO TO 20
      IF (MOD(I(2),3).NE.0) GO TO 20
      IF (MOD(I(3),3).NE.0) GO TO 20
      K=1
      RETURN
   20 IF (IMAX.LT.5) RETURN
      IF (MOD(I(1),5).NE.0) GO TO 30
      IF (MOD(I(2),5).NE.0) GO TO 30
      IF (MOD(I(3),5).NE.0) GO TO 30
      K=1
      RETURN
   30 IF (IMAX.LT.7) RETURN
      IF (MOD(I(1),7).NE.0) GO TO 40
      IF (MOD(I(2),7).NE.0) GO TO 40
      IF (MOD(I(3),7).NE.0) GO TO 40
      K=1
      RETURN
   40 IF (IMAX.LT.11) RETURN
      IF (MOD(I(1),11).NE.0) GO TO 50
      IF (MOD(I(2),11).NE.0) GO TO 50
      IF (MOD(I(3),11).NE.0) GO TO 50
      K=1
      RETURN
   50 IF (IMAX.LT.13) RETURN
      IF (MOD(I(1),13).NE.0) GO TO 60
      IF (MOD(I(2),13).NE.0) GO TO 60
      IF (MOD(I(3),13).NE.0) GO TO 60
      K=1
      RETURN
   60 IF (IMAX.LT.17) RETURN
      IF (MOD(I(1),17).NE.0) GO TO 70
      IF (MOD(I(2),17).NE.0) GO TO 70
      IF (MOD(I(3),17).NE.0) GO TO 70
      K=1
      RETURN
   70 IF (IMAX.LT.19) RETURN
      IF (MOD(I(1),19).NE.0) RETURN
      IF (MOD(I(2),19).NE.0) RETURN
      IF (MOD(I(3),19).NE.0) RETURN
      K=1
      RETURN
      END
      REAL*8 FUNCTION DSIND(X)
      REAL*8 X,DEGPI,DEG,PI
      DATA DEG,PI/180.D00,3.14159265358979D00/
      DEGPI=DEG/PI
      DSIND=DSIN(X/DEGPI)
      RETURN
      END
      REAL*8 FUNCTION DCOSD(X)
      REAL*8 X,DEGPI,DEG,PI
      DATA DEG,PI/180.D00,3.14159265358979D00/
      DEGPI=DEG/PI
      DCOSD=DCOS(X/DEGPI)
      RETURN
      END
      SUBROUTINE DOTP (A1, A2, A3, B1, B2, B3, ANG)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X,DEGPI,DEG,PI
      DATA DEG,PI/180.D00,3.14159265358979D00/
      DEGPI=DEG/PI
      ANORM=SQRT(A1**2+A2**2+A3**2)
      BNORM=SQRT(B1**2+B2**2+B3**2)
      AN1=A1/ANORM
      AN2=A2/ANORM
      AN3=A3/ANORM
      BN1=B1/BNORM
      BN2=B2/BNORM
      BN3=B3/BNORM
      CANG=AN1*BN1+AN2*BN2+AN3*BN3
C     SOLVE ROUNDOFF ERROR OCCASIONAL ISSUE
      IF (cang.gt.1.000000.and.cang.le.1.0000001) cang=1.0
      IF (cang.lt.-1.000000.and.cang.ge.-1.0000001) cang=-1.0
      ANG=DEGPI*DACOS(CANG)
      RETURN
      END
      SUBROUTINE SORT (JTEST, MULT, ATMOD, NEWOLD, ITEST)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /NEAT/ JD(3,500,64),JR(3,500,64),JDT(3,500,64),JRT(3,500,
     164),JMULT(500,64),JTWIN(500,64),OBL(500,64),JC(64)
      COMMON /FINAL/ JDN(3,500,32),JRN(3,500,32),JDTN(3,500,32),JRTN(3,
     1500,32),JMULTN(500,32),JTWINN(500,32),OBLN(500,32),JCNEW(32),
     2NTEST(500),NLAWS
      DIMENSION INDD(3), INDR(3)
      CHARACTER ATMOD*1
C 
      IF (JTEST.EQ.0) RETURN
      IF (NEWOLD.EQ.2) GO TO 40
      DO 30 K=1,jtest-1
        DO 20 L=K+1,jtest
          IF (OBL(K,MULT).GT.OBL(L,MULT)) THEN
            TEMP=OBL(K,MULT)
            OBL(K,MULT)=OBL(L,MULT)
            OBL(L,MULT)=TEMP
            ITEMP=JTWIN(K,MULT)
            JTWIN(K,MULT)=JTWIN(L,MULT)
            JTWIN(L,MULT)=ITEMP
            ITEMP=JMULT(K,MULT)
            JMULT(K,MULT)=JMULT(L,MULT)
            JMULT(L,MULT)=ITEMP
            DO 10 II=1,3
              INDD(II)=JD(II,K,MULT)
              JD(II,K,MULT)=JD(II,L,MULT)
              JD(II,L,MULT)=INDD(II)
              INDR(II)=JR(II,K,MULT)
              JR(II,K,MULT)=JR(II,L,MULT)
              JR(II,L,MULT)=INDR(II)
C              IF(ATMOD.EQ.'P'.OR.ATMOD.EQ.'p') GO TO 12
              IF (ITEST.EQ.0) GO TO 10
              INDD(II)=JDT(II,K,MULT)
              JDT(II,K,MULT)=JDT(II,L,MULT)
              JDT(II,L,MULT)=INDD(II)
              INDR(II)=JRT(II,K,MULT)
              JRT(II,K,MULT)=JRT(II,L,MULT)
              JRT(II,L,MULT)=INDR(II)
   10       CONTINUE
          END IF
   20   CONTINUE
   30 CONTINUE
      GO TO 80
   40 DO 70 K=1,jtest-1
        DO 60 L=K+1,jtest
          IF (OBLN(K,MULT).GT.OBLN(L,MULT)) THEN
            TEMP=OBLN(K,MULT)
            OBLN(K,MULT)=OBLN(L,MULT)
            OBLN(L,MULT)=TEMP
            ITEMP=JTWINN(K,MULT)
            JTWINN(K,MULT)=JTWINN(L,MULT)
            JTWINN(L,MULT)=ITEMP
            ITEMP=JMULTN(K,MULT)
            JMULTN(K,MULT)=JMULTN(L,MULT)
            JMULTN(L,MULT)=ITEMP
            DO 50 II=1,3
              INDD(II)=JDN(II,K,MULT)
              JDN(II,K,MULT)=JDN(II,L,MULT)
              JDN(II,L,MULT)=INDD(II)
              INDR(II)=JRN(II,K,MULT)
              JRN(II,K,MULT)=JRN(II,L,MULT)
              JRN(II,L,MULT)=INDR(II)
C              IF(ATMOD.EQ.'P'.OR.ATMOD.EQ.'p') GO TO 32
              IF (ITEST.EQ.0) GO TO 50
              INDD(II)=JDTN(II,K,MULT)
              JDTN(II,K,MULT)=JDTN(II,L,MULT)
              JDTN(II,L,MULT)=INDD(II)
              INDR(II)=JRTN(II,K,MULT)
              JRTN(II,K,MULT)=JRTN(II,L,MULT)
              JRTN(II,L,MULT)=INDR(II)
   50       CONTINUE
          END IF
   60   CONTINUE
   70 CONTINUE
   80 RETURN
      END
      SUBROUTINE COPRIME2 (AIHXT, AIPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AIHXT(3), AIPT(3)
      DATA Z2/2.0D00/
      REAL*8 IFAC
   10 IFAC=1.
      IF (DMOD(AIHXT(1),Z2).NE.0..OR.DMOD(AIHXT(2),Z2)
     1.NE.0..OR.DMOD(AIHXT(3),Z2).NE.0.) GO TO 30
      IFAC=2.
      DO 20 I=1,3
   20   AIHXT(I)=AIHXT(I)/IFAC
      GO TO 10
C 
   30 IFAC=1.
      IF (DMOD(AIPT(1),Z2).NE.0..OR.DMOD(AIPT(2),Z2)
     1.NE.0..OR.DMOD(AIPT(3),Z2).NE.0.) GO TO 50
      IFAC=2.
      DO 40 I=1,3
   40   AIPT(I)=AIPT(I)/IFAC
      GO TO 30
   50 RETURN
      END
      SUBROUTINE REGROUP (IMAX, ATMOD, IOUT, ITEST)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /NEAT/ JD(3,500,64),JR(3,500,64),JDT(3,500,64),JRT(3,500,
     164),JMULT(500,64),JTWIN(500,64),OBL(500,64),JC(64)
      COMMON /FINAL/ JDN(3,500,32),JRN(3,500,32),JDTN(3,500,32),JRTN(3,
     1500,32),JMULTN(500,32),JTWINN(500,32),OBLN(500,32),JCNEW(32),
     2NTEST(500),NLAWS
      COMMON /SYST/ XSYS
      CHARACTER ATMOD*1,XSYS*4
C 
C       JC COUNTS NUMBER IN A GIVEN DOT, JCNEW COUNTS NUMBER IN A GIVEN
C 
C       GET TOTALS FOR N ODD (DOT) AND 2N (TWIN INDEX N); JJ OVER ALL TW
      IF (MOD(IMAX,2).NE.0) IMAX=IMAX-1
      DO 10 JJ=1,IMAX/2,2
   10   JCNEW(JJ)=JC(JJ)+JC(JJ*2)
      DO 20 JJ=2,IMAX/2,2
   20   JCNEW(JJ)=JC(2*JJ)
      DO 70 I=1,IMAX/2,2
        JJ=I
        DO 40 K=1,JC(JJ)
          IF (JC(I).EQ.0) GO TO 40
          JMULTN(K,I)=JMULT(K,I)
          JTWINN(K,I)=JTWIN(K,I)
          DO 30 L=1,3
            JDN(L,K,I)=JD(L,K,I)
            JRN(L,K,I)=JR(L,K,I)
            JDTN(L,K,I)=JDT(L,K,I)
   30       JRTN(L,K,I)=JRT(L,K,I)
          OBLN(K,I)=OBL(K,I)
   40   CONTINUE
        DO 60 K=JC(JJ)+1,K+JC(JJ*2)
          JMULTN(K,I)=JMULT(K-JC(JJ),2*I)
          JTWINN(K,I)=JTWIN(K-JC(JJ),2*I)
          DO 50 L=1,3
            JDN(L,K,I)=JD(L,K-JC(JJ),2*I)
            JRN(L,K,I)=JR(L,K-JC(JJ),2*I)
            JDTN(L,K,I)=JDT(L,K-JC(JJ),2*I)
   50       JRTN(L,K,I)=JRT(L,K-JC(JJ),2*I)
   60     OBLN(K,I)=OBL(K-JC(JJ),2*I)
   70 CONTINUE
      DO 90 I=2,IMAX/2,2
        DO 90 K=1,JCNEW(I)
        IF (JCNEW(I).EQ.0) GO TO 90
        JMULTN(K,I)=JMULT(K,2*I)
        JTWINN(K,I)=JTWIN(K,2*I)
        DO 80 L=1,3
          JDN(L,K,I)=JD(L,K,2*I)
          JRN(L,K,I)=JR(L,K,2*I)
          JDTN(L,K,I)=JDT(L,K,2*I)
   80     JRTN(L,K,I)=JRT(L,K,2*I)
        OBLN(K,I)=OBL(K,2*I)
   90 CONTINUE
      IF (NVERBOSE.EQ.1) GO TO 160
      WRITE (IOUT,100)
      WRITE (*,100)
  100 FORMAT (//5X,'OBLIQUE-STYLE OUTPUT, SORTED ON ''M'' ','(TWIN INDEX
     1) FOLLOWS: '/)
C       IF(ATMOD.NE.'P'.AND.ATMOD.NE.'p') GO TO 43
      IF (ITEST.EQ.1) GO TO 130
      WRITE (IOUT,110)
      WRITE (*,110)
  110 FORMAT (/15X,'ROWS',22X,'PRODUCTS'//11X,'BUERGER CELL')
      WRITE (IOUT,120)
      WRITE (*,120)
  120 FORMAT (/,7X,'DIRECT',7X,'RECIPROCAL',9X,'DOT',4X,'M',4X,' OMEGA')
      GO TO 160
  130 WRITE (IOUT,140)
      WRITE (*,140)
  140 FORMAT (/15X,'ROWS',22X,'PRODUCTS'//11X,'BUERGER CELL',42X,'      
     1INPUT CELL')
      WRITE (IOUT,150)
      WRITE (*,150)
  150 FORMAT (/,7X,'DIRECT',7X,'RECIPROCAL',9X,'DOT',4X,'M',4X,' OMEGA',
     18X,'DIRECT',8X,'RECIPROCAL')
  160 DO 190 J=1,IMAX/2
        NEWOLD=2
        CALL SORT (JCNEW(J), J, ATMOD, NEWOLD, ITEST)
        DO 180 K=1,JCNEW(J)
          IF (JCNEW(J).EQ.0) GO TO 180
C       IF(ATMOD.NE.'P'.AND.ATMOD.NE.'p') GO TO 22
          IF (ITEST.EQ.1) GO TO 170
          IF (NVERBOSE.EQ.1) GO TO 180
          WRITE (IOUT,370) (JDN(L,K,J),L=1,3),(JRN(L,K,J),L=1,3),
     1     JMULTN(K,J),JTWINN(K,J),OBLN(K,J)
          WRITE (*,370) (JDN(L,K,J),L=1,3),(JRN(L,K,J),L=1,3),JMULTN(K,
     1     J),JTWINN(K,J),OBLN(K,J)
          GO TO 180
C       JDN AND JRN = BUERGER CELL DIR & RECIP TWIN DIRECTIONS, RESP.
C 
C       JDTN AND JRTN = INPUT CELL DIR & RECIP TWIN DIRECIONS, RESP.
C 
C       JMULT = SCALAR ; JTWINN = TWIN INDEX ; OBLN = OBLIQUITY
C 
  170     WRITE (IOUT,380) (JDN(L,K,J),L=1,3),(JRN(L,K,J),L=1,3),
     1     JMULTN(K,J),JTWINN(K,J),OBLN(K,J),(JDTN(L,K,J),L=1,3),
     2     (JRTN(L,K,J),L=1,3)
          WRITE (*,380) (JDN(L,K,J),L=1,3),(JRN(L,K,J),L=1,3),JMULTN(K,
     1     J),JTWINN(K,J),OBLN(K,J),(JDTN(L,K,J),L=1,3),(JRTN(L,K,J),L=
     2     1,3)
  180   CONTINUE
  190 CONTINUE
  
      IF (XSYS.EQ.'    ') RETURN
      IF (ITEST.EQ.0) GO TO 280
      IF (XSYS.EQ.'MONO') GO TO 200
      WRITE (*,390)
      WRITE (IOUT,390)
      WRITE (*,140)
      WRITE (*,150)
      WRITE (IOUT,140)
      WRITE (IOUT,150)
      GO TO 210
  200 WRITE (*,400)
      WRITE (IOUT,400)
      WRITE (*,140)
      WRITE (*,150)
      WRITE (IOUT,140)
      WRITE (IOUT,150)
  210 DO 270 J=1,IMAX/2
        DO 220 I=1,500
  220     NTEST(I)=0
        DO 260 K=1,JCNEW(J)
          IF (NTEST(K).EQ.0) WRITE (IOUT,380) (JDN(L,K,J),L=1,3),(JRN(L,
     1     K,J),L=1,3),JMULTN(K,J),JTWINN(K,J),OBLN(K,J),(JDTN(L,K,J),L=
     2     1,3),(JRTN(L,K,J),L=1,3)
          IF (NTEST(K).EQ.0) WRITE (*,380) (JDN(L,K,J),L=1,3),(JRN(L,K,
     1     J),L=1,3),JMULTN(K,J),JTWINN(K,J),OBLN(K,J),(JDTN(L,K,J),L=1,
     2     3),(JRTN(L,K,J),L=1,3)
          IF (NTEST(K).EQ.0) CALL TLCALC (JRTN(1,K,J), JRTN(2,K,J),
     1     JRTN(3,K,J), IOUT)
          IF (NTEST(K).EQ.0) CALL TLCALC2 (JDTN(1,K,J), JDTN(2,K,J),
     1     JDTN(3,K,J), IOUT2)
          IF (NTEST(K).EQ.1) GO TO 260
          IF (XSYS.EQ.'MONO') GO TO 240
          DO 230 JJ=1,3
            IF (JMULTN(K,J).NE.JMULTN(K+JJ,J)) GO TO 260
            IF (JDTN(1,K,J).EQ.-JDTN(1,K+1,J).AND.JDTN(2,K,J).EQ.-
     1       JDTN(2,K+JJ,J).AND.JDTN(3,K,J).EQ.-JDTN(3,K+JJ,J)
     2       .AND.JRTN(1,K,J).EQ.-JRTN(1,K+JJ,J).AND.JRTN(2,K,J).EQ.-
     3       JRTN(2,K+JJ,J).AND.JRTN(3,K,J).EQ.-JRTN(3,K+JJ,J)) NTEST(K+
     4       JJ)=1
  230     CONTINUE
          GO TO 260
  240     DO 250 JJ=1,3
            IF (JMULTN(K,J).NE.JMULTN(K+JJ,J)) GO TO 260
            IF (JDTN(1,K,J).EQ.-JDTN(1,K+JJ,J).AND.JDTN(2,K,J)
     1       .EQ.JDTN(2,K+JJ,J).AND.JDTN(3,K,J).EQ.-JDTN(3,K+JJ,J)
     2       .AND.JRTN(1,K,J).EQ.-JRTN(1,K+JJ,J).AND.JRTN(2,K,J)
     3       .EQ.JRTN(2,K+JJ,J).AND.JRTN(3,K,J).EQ.-JRTN(3,K+JJ,J))
     4       NTEST(K+JJ)=1
            IF (JDTN(1,K,J).EQ.JDTN(1,K+JJ,J).AND.JDTN(2,K,J).EQ.-
     1       JDTN(2,K+JJ,J).AND.JDTN(3,K,J).EQ.JDTN(3,K+JJ,J)
     2       .AND.JRTN(1,K,J).EQ.JRTN(1,K+JJ,J).AND.JRTN(2,K,J).EQ.-
     3       JRTN(2,K+JJ,J).AND.JRTN(3,K,J).EQ.JRTN(3,K+JJ,J)) NTEST(K+
     4       JJ)=1
            IF (JDTN(1,K,J).EQ.-JDTN(1,K+1,J).AND.JDTN(2,K,J).EQ.-
     1       JDTN(2,K+JJ,J).AND.JDTN(3,K,J).EQ.-JDTN(3,K+JJ,J)
     2       .AND.JRTN(1,K,J).EQ.-JRTN(1,K+JJ,J).AND.JRTN(2,K,J).EQ.-
     3       JRTN(2,K+JJ,J).AND.JRTN(3,K,J).EQ.-JRTN(3,K+JJ,J)) NTEST(K+
     4       JJ)=1
  250     CONTINUE
  260   CONTINUE
  270 CONTINUE
      CALL TLCALC (0, 0, 0, IOUT)
      CALL TLCALC2 (0, 0, 0, IOUT)
      RETURN
  280 IF (XSYS.EQ.'MONO') GO TO 290
      WRITE (*,390)
      WRITE (IOUT,390)
      WRITE (*,110)
      WRITE (*,120)
      WRITE (IOUT,110)
      WRITE (IOUT,120)
      GO TO 300
  290 WRITE (*,400)
      WRITE (IOUT,400)
      WRITE (*,110)
      WRITE (*,120)
      WRITE (IOUT,110)
      WRITE (IOUT,120)
  300 DO 360 J=1,IMAX/2
        DO 310 I=1,500
  310     NTEST(I)=0
        DO 350 K=1,JCNEW(J)
          IF (NTEST(K).EQ.0) WRITE (IOUT,370) (JDN(L,K,J),L=1,3),(JRN(L,
     1     K,J),L=1,3),JMULTN(K,J),JTWINN(K,J),OBLN(K,J)
          IF (NTEST(K).EQ.0) WRITE (*,370) (JDN(L,K,J),L=1,3),(JRN(L,K,
     1     J),L=1,3),JMULTN(K,J),JTWINN(K,J),OBLN(K,J)
          IF (NTEST(K).EQ.0) CALL TLCALC (JRN(1,K,J), JRN(2,K,J), JRN(3,
     1     K,J), IOUT)
          IF (NTEST(K).EQ.0) CALL TLCALC2 (JDN(1,K,J), JDN(2,K,J),
     1     JDN(3,K,J), IOUT2)
          IF (NTEST(K).EQ.1) GO TO 350
          IF (XSYS.EQ.'MONO') GO TO 330
          DO 320 JJ=1,3
            IF (JMULTN(K,J).NE.JMULTN(K+JJ,J)) GO TO 350
            IF (JDN(1,K,J).EQ.-JDN(1,K+1,J).AND.JDN(2,K,J).EQ.-JDN(2,K+
     1       JJ,J).AND.JDN(3,K,J).EQ.-JDN(3,K+JJ,J).AND.JRN(1,K,J).EQ.-
     2       JRN(1,K+JJ,J).AND.JRN(2,K,J).EQ.-JRN(2,K+JJ,J).AND.JRN(3,K,
     3       J).EQ.-JRN(3,K+JJ,J)) NTEST(K+JJ)=1
  320     CONTINUE
          GO TO 350
  330     DO 340 JJ=1,3
            IF (JMULTN(K,J).NE.JMULTN(K+JJ,J)) GO TO 350
            IF (JDN(1,K,J).EQ.-JDN(1,K+JJ,J).AND.JDN(2,K,J).EQ.JDN(2,K+
     1       JJ,J).AND.JDN(3,K,J).EQ.-JDN(3,K+JJ,J).AND.JRN(1,K,J).EQ.-
     2       JRN(1,K+JJ,J).AND.JRN(2,K,J).EQ.JRN(2,K+JJ,J).AND.JRN(3,K,
     3       J).EQ.-JRN(3,K+JJ,J)) NTEST(K+JJ)=1
            IF (JDN(1,K,J).EQ.JDN(1,K+JJ,J).AND.JDN(2,K,J).EQ.-JDN(2,K+
     1       JJ,J).AND.JDN(3,K,J).EQ.JDN(3,K+JJ,J).AND.JRN(1,K,J)
     2       .EQ.JRN(1,K+JJ,J).AND.JRN(2,K,J).EQ.-JRN(2,K+JJ,J)
     3       .AND.JRN(3,K,J).EQ.JRN(3,K+JJ,J)) NTEST(K+JJ)=1
            IF (JDN(1,K,J).EQ.-JDN(1,K+1,J).AND.JDN(2,K,J).EQ.-JDN(2,K+
     1       JJ,J).AND.JDN(3,K,J).EQ.-JDN(3,K+JJ,J).AND.JRN(1,K,J).EQ.-
     2       JRN(1,K+JJ,J).AND.JRN(2,K,J).EQ.-JRN(2,K+JJ,J).AND.JRN(3,K,
     3       J).EQ.-JRN(3,K+JJ,J)) NTEST(K+JJ)=1
  340     CONTINUE
  350   CONTINUE
  360 CONTINUE
  
  
  370 FORMAT (2X,'[',3I4,']',2X,'(',3I4,')',I10,I5,F10.4)
  380 FORMAT (2X,'[',3I4,']',2X,'(',3I4,')',I10,I5,F10.4,4X,'[',3I4,']',
     12x,'(',3I4,')')
  390 FORMAT (///,'  SYMMETRY-INDEPENDENT DIRECT AND RECIPROCAL SETS',' 
     1FOR TRICLINIC SYSTEM FOLLOW'//)
  400 FORMAT (///,'  SYMMETRY-INDEPENDENT DIRECT AND RECIPROCAL SETS',' 
     1FOR MONOCLINIC SYSTEM FOLLOW'//)
      CALL TLCALC (0, 0, 0, IOUT)
      CALL TLCALC2 (0, 0, 0, IOUT)
      RETURN
      END
      SUBROUTINE STANDCALC (CELL, CELLA, IOUT)
C 
C     CALCULATE METRIC TENSOR, ROLLETT L AND U MATRICES FOR TWIN LAW CAL
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 HR,KR,LR,L,LT
      COMMON /STD/ L(3,3),LT(3,3),U(3,3),G(3,3),GINV(3,3),UT(3,3)
      DIMENSION CELL(3), CELLA(3)
      COMMON /NEAT/ JD(3,500,64),JR(3,500,64),JDT(3,500,64),JRT(3,500,
     164),JMULT(500,64),JTWIN(500,64),OBL(500,64),JC(64)
      COMMON /FINAL/ JDN(3,500,32),JRN(3,500,32),JDTN(3,500,32),JRTN(3,
     1500,32),JMULTN(500,32),JTWINN(500,32),OBLN(500,32),JCNEW(32),
     2NTEST(500),NLAWS
      COMMON /SYST/ XSYS
      CHARACTER ATMOD*1,XSYS*4
      DATA DTOA/57.29578/
C 
C 
C     ORTHOGONALISATION AS DONE IN ROLLETT, COMPUTING METHODS IN
C     CRYSTALLOGRAPHY, PERGAMON PRESS, 1965, PP. 22-24.
C     WRITE (*,500)
C     READ(*,*) A, B, C, ALP, BET, GAM
C     OPEN(12,STATUS='UNKNOWN',FILE='CTOOL.OUT',ACCESS='SEQUENTIAL',
C    1FORM='FORMATTED')
      A=CELL(1)
      B=CELL(2)
      C=CELL(3)
      ALP=CELLA(1)
      BET=CELLA(2)
      GAM=CELLA(3)
      WRITE (IOUT,10) A,B,C,ALP,BET,GAM
      ALPR=ALP/DTOA
      BETR=BET/DTOA
      GAMR=GAM/DTOA
      COSA=COS(ALPR)
      COSB=COS(BETR)
      COSG=COS(GAMR)
      SINA=SQRT(1.-COSA**2)
      SINB=SQRT(1.-COSB**2)
      SING=SQRT(1.-COSG**2)
      COSGST=(COSA*COSB-COSG)/(SINA*SINB)
      SINGST=SQRT(1.-COSGST**2)
      VOL=A*B*C*SQRT(1.-COSA**2-COSB**2-COSG**2+2.*COSA*COSB*COSG)
      ASTAR=B*C*SINA/VOL
      BSTAR=A*C*SINB/VOL
      CSTAR=A*B*SING/VOL
      COSAST=(COSB-COSG-COSA)/(SINB*SING)
      COSBST=(COSA*COSG-COSB)/(SINA*SING)
      ALSTAR=DTOA*ACOS(COSAST)
      BESTAR=DTOA*ACOS(COSBST)
      GASTAR=DTOA*ACOS(COSGST)
      WRITE (IOUT,20) ASTAR,BSTAR,CSTAR,ALSTAR,BESTAR,GASTAR
      L(1,1)=A*SINB*SINGST
      L(1,2)=0.
      L(1,3)=0.
      L(2,1)=-1.*A*SINB*COSGST
      L(2,2)=B*SINA
      L(2,3)=0.
      L(3,1)=A*COSB
      L(3,2)=B*COSA
      L(3,3)=C
      WRITE (IOUT,30) ((L(I,J),J=1,3),I=1,3)
C     TRANSPOSE OF L MATRIX
      CALL MATRANS (L, LT)
C     INVERT L^T
      CALL MATINV (LT, U, D)
      CALL MATRANS (U, UT)
      WRITE (IOUT,40) ((U(I,J),J=1,3),I=1,3)
C     SET UP METRIC TENSOR AND ITS INVERSE
      G(1,1)=A**2
      G(1,2)=A*B*COSG
      G(1,3)=A*C*COSB
      G(2,1)=G(1,2)
      G(2,2)=B**2
      G(2,3)=B*C*COSA
      G(3,1)=G(1,3)
      G(3,2)=G(2,3)
      G(3,3)=C**2
      WRITE (IOUT,50) ((G(I,J),J=1,3),I=1,3)
C     GET INVERSE OF METRIC TENSOR
      CALL MATINV (G, GINV, D)
   10 FORMAT (//' DIRECT CELL CONSTANTS: ',6F10.4)
   20 FORMAT (//' RECIPROCAL CELL CONSTANTS:',3F10.6,3F10.4)
   30 FORMAT (//' DIRECT ORTHOGONALISATION MATRIX, J. S. ROLLETT,'' COMP
     1UTING METHODS IN CRYSTALLOGRAPHY,'/' PERGAMON, 1965, PP. ','22-24:
     2 ',//3(3F12.4,/))
   40 FORMAT (//' RECIP. ORTHOGONALISATION MATRIX, J. S. ROLLETT,'' COMP
     1UTING METHODS IN CRYSTALLOGRAPHY,'/' PERGAMON, 1965, PP. ','22-24:
     2 ',//3(3F12.6,/))
   50 FORMAT (//' METRIC TENSOR',//3(3F12.6,/))
      RETURN
      END
      SUBROUTINE TLCALC (IH, IK, IL, IOUT)
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 HR,KR,LR,L,LT,HDSV,KDSV,LDSV,HRSV,KRSV,LRSV,DSV,D,DP,RSV,
     1RP,R
      COMMON /STD/ L(3,3),LT(3,3),U(3,3),G(3,3),GINV(3,3),UT(3,3)
      REAL*8 CELL(3),CELLA(3),TEMP(3,3)
      COMMON /NEAT/ JD(3,500,64),JR(3,500,64),JDT(3,500,64),JRT(3,500,
     164),JMULT(500,64),JTWIN(500,64),OBL(500,64),JC(64)
      COMMON /FINAL/ JDN(3,500,32),JRN(3,500,32),JDTN(3,500,32),JRTN(3,
     1500,32),JMULTN(500,32),JTWINN(500,32),OBLN(500,32),JCNEW(32),
     2NTEST(500),NLAWS,NLAWS2
      COMMON /TLAW/ HDSV(1500),KDSV(1500),LDSV(1500),DP(3,3),D(3,3),
     1HRSV(1500),KRSV(1500),LRSV(1500),DSV(3,3,1500),RP(3,3),R(3,3),
     2RSV(3,3,1500)
      COMMON /SYST/ XSYS
      CHARACTER ATMOD*1,XSYS*4
C 
C     READ(*,*) HR(1),KR(1),LR(1)
C     IF(HR(1).EQ.100) GO TO 220
      IF (IH.EQ.0.AND.IK.EQ.0.AND.IL.EQ.0) GO TO 20
      NLAWS=NLAWS+1
      NLAWS2=NLAWS
      HR=IH
      KR=IK
      LR=IL
C     ORTHOGONALISE (hkl)
      C1R=U(1,1)*HR+U(1,2)*KR+U(1,3)*LR
      C2R=U(2,1)*HR+U(2,2)*KR+U(2,3)*LR
      C3R=U(3,1)*HR+U(3,2)*KR+U(3,3)*LR
      CRL=SQRT(C1R**2+C2R**2+C3R**2)
      VV1=C1R/CRL
      VV2=C2R/CRL
      VV3=C3R/CRL
      RP(1,1)=2.*VV1**2-1.
      RP(2,2)=2.*VV2**2-1.
      RP(3,3)=2.*vv3**2-1.
      RP(1,2)=2.*VV1*VV2
      RP(1,3)=2.*VV1*VV3
      RP(2,3)=2.*VV2*VV3
      RP(2,1)=RP(1,2)
      RP(3,1)=RP(1,3)
      RP(3,2)=RP(2,3)
C     GENERATE TWIN LAW MATRIX FOR (hkl)
      CALL MATMULT (RP, U, 3, TEMP)
      CALL MATMULT (LT, TEMP, 3, R)
      HRSV(NLAWS)=HR
      KRSV(NLAWS)=KR
      LRSV(NLAWS)=LR
      DO 10 I=1,3
        DO 10 J=1,3
   10   RSV(I,J,NLAWS)=R(I,J)
      GO TO 60
   20 DO 30 JJ=1,NLAWS
        WRITE (*,40) HRSV(JJ),KRSV(JJ),LRSV(JJ)
        WRITE (IOUT,40) HRSV(JJ),KRSV(JJ),LRSV(JJ)
        WRITE (*,50) ((RSV(I,J,JJ),J=1,3),I=1,3)
   30   WRITE (IOUT,50) ((RSV(I,J,JJ),J=1,3),I=1,3)
   40 FORMAT (//' TWIN LAW FOR 180 DEG ROTATION ABOUT R.L. VECTOR: ',
     13F8.3)
   50 FORMAT (//,3(3F10.6,/))
   60 RETURN
      END
      SUBROUTINE TLCALC2 (IH, IK, IL, IOUT)
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 HR,KR,LR,L,LT,HDSV,KDSV,LDSV,HRSV,KRSV,LRSV,DSV,D,DP,RSV,
     1RP,R
      COMMON /STD/ L(3,3),LT(3,3),U(3,3),G(3,3),GINV(3,3),UT(3,3)
      REAL*8 CELL(3),CELLA(3),TEMP(3,3)
      COMMON /NEAT/ JD(3,500,64),JR(3,500,64),JDT(3,500,64),JRT(3,500,
     164),JMULT(500,64),JTWIN(500,64),OBL(500,64),JC(64)
      COMMON /FINAL/ JDN(3,500,32),JRN(3,500,32),JDTN(3,500,32),JRTN(3,
     1500,32),JMULTN(500,32),JTWINN(500,32),OBLN(500,32),JCNEW(32),
     2NTEST(500),NLAWS,NLAWS2
      COMMON /TLAW/ HDSV(1500),KDSV(1500),LDSV(1500),DP(3,3),D(3,3),
     1HRSV(1500),KRSV(1500),LRSV(1500),DSV(3,3,1500),RP(3,3),R(3,3),
     2RSV(3,3,1500)
      COMMON /SYST/ XSYS
      CHARACTER ATMOD*1,XSYS*4
C 
C     READ(*,*) HR(1),KR(1),LR(1)
C     IF(HR(1).EQ.100) GO TO 220
      IF (IH.EQ.0.AND.IK.EQ.0.AND.IL.EQ.0) GO TO 20
      NLAWS=NLAWS2
      HD=IH
      KD=IK
      LD=IL
C     ORTHOGONALISE [hkl]
      C1D=L(1,1)*HD+L(1,2)*KD+L(1,3)*LD
      C2D=L(2,1)*HD+L(2,2)*KD+L(2,3)*LD
      C3D=L(3,1)*HD+L(3,2)*KD+L(3,3)*LD
      CDL=SQRT(C1D**2+C2D**2+C3D**2)
      VV1=C1D/CDL
      VV2=C2D/CDL
      VV3=C3D/CDL
      DP(1,1)=2.*VV1**2-1.
      DP(2,2)=2.*VV2**2-1.
      DP(3,3)=2.*vv3**2-1.
      DP(1,2)=2.*VV1*VV2
      DP(1,3)=2.*VV1*VV3
      DP(2,3)=2.*VV2*VV3
      DP(2,1)=DP(1,2)
      DP(3,1)=DP(1,3)
      DP(3,2)=DP(2,3)
C     GENERATE TWIN LAW MATRIX FOR [hkl]
      CALL MATMULT (DP, U, 3, TEMP)
      CALL MATMULT (LT, TEMP, 3, D)
      HDSV(NLAWS)=HD
      KDSV(NLAWS)=KD
      LDSV(NLAWS)=LD
      DO 10 I=1,3
        DO 10 J=1,3
   10   DSV(I,J,NLAWS)=D(I,J)
      GO TO 60
   20 DO 30 JJ=1,NLAWS
        WRITE (*,40) HDSV(JJ),KDSV(JJ),LDSV(JJ)
        WRITE (IOUT,40) HDSV(JJ),KDSV(JJ),LDSV(JJ)
        WRITE (*,50) ((DSV(I,J,JJ),J=1,3),I=1,3)
   30   WRITE (IOUT,50) ((DSV(I,J,JJ),J=1,3),I=1,3)
   40 FORMAT (//' TWIN LAW FOR 180 DEG ROTATION ABOUT DIRECT VECTOR: ',
     13F8.3)
   50 FORMAT (//,3(3F10.6,/))
   60 RETURN
      END
      SUBROUTINE MATRANS (A, AT)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(3,3),AT(3,3)
      DO 20 I=1,3
        DO 10 J=1,3
          AT(J,I)=A(I,J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE MATINV (A, B, D)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER*2 K,N
c      REAL*4    A,B,D
      DIMENSION A(3,3), B(3,3)
      CALL VMULT (A(1,2), A(1,3), B(1,1))
      CALL VMULT (A(1,3), A(1,1), B(1,2))
      CALL VMULT (A(1,1), A(1,2), B(1,3))
      D=A(1,1)*B(1,1)+A(2,1)*B(2,1)+A(3,1)*B(3,1)
      IF (D.EQ.0.0) RETURN
      N=1
   10 K=1
   20 B(K,N)=B(K,N)/D
      K=K+1
      IF (K.LE.3) GO TO 20
      N=N+1
      IF (N.LE.3) GO TO 10
      CALL MACOL (B)
      RETURN
      END
      SUBROUTINE matmult (a, b, n, c)
      REAL*8 a(n,n),b(n,n),c(n,n)
      DO 10 i=1,n
        DO 10 j=1,n
        c(i,j)=0.
        DO 10 k=1,n
   10   c(i,j)=c(i,j)+a(i,k)*b(k,j)
      RETURN
      END
      SUBROUTINE VMULT (A, B, C)
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION A(3), B(3), C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
      SUBROUTINE MACOL (A)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER*2 K,N
C     REAL*4    A,T
      DIMENSION A(3,3)
      N=1
   10 K=2
   20 T=A(K,N)
      A(K,N)=A(N,K)
      A(N,K)=T
      IF (N.EQ.3) RETURN
      K=K+1
      IF (K.LE.3) GO TO 20
      N=3
      GO TO 10
      END
