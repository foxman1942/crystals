      PROGRAM CELLTOOL
      implicit real*8 (a-h,o-z)
      REAL*8 L(3,3),LT(3,3),U(3,3),HD(50),KD(50),LD(50),
     1HR(50),KR(50),LR(50),C1D(50),C2D(50),C3D(50),
     2C1R(50),C2R(50),C3R(50),ANG(50),HR2(50),KR2(50),LR2(50),
     3C1R2(50),C2R2(50),C3R2(50),HD2(50),KD2(50),LD2(50),
     4C1D2(50),C2D2(50),C3D2(50),LEN(50),P1(50),P2(50),P3(50),
     5NN1(50),NN2(50),NN3(50),VV1(50),VV2(50),VV3(50),
     6UU1(50),UU2(50),UU3(50),G(3,3),GINV(3,3),R(3,3),RP(3,3),
     7UVW(3),UVWP(3),UVWPM(3),HKL(3),HKLP(3),HKLPM(3),TEMP(3,3),
     8VVD(3,3),RPD(3,3),RD(3,3),VV1D(50),VV2D(50),VV3D(50),UT(3,3)
      DATA DTOA/57.29578/
      CHARACTER*4 TEST
C     ORTHOGONALISATION AS DONE IN ROLLETT, COMPUTING METHODS IN 
C     CRYSTALLOGRAPHY, PERGAMON PRESS, 1965, PP. 22-24.
      WRITE (*,500)
      READ(*,*) A, B, C, ALP, BET, GAM
      OPEN(12,STATUS='UNKNOWN',FILE='CTOOL.OUT',ACCESS='SEQUENTIAL',
     1FORM='FORMATTED')
      write(12,501) A,B,C,ALP,BET,GAM
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
      COSAST=(COSB*COSG-COSA)/(SINB*SING)
      COSBST=(COSA*COSG-COSB)/(SINA*SING)
      ALSTAR=DTOA*ACOS(COSAST)
      BESTAR=DTOA*ACOS(COSBST)
      GASTAR=DTOA*ACOS(COSGST)
      WRITE(12,502) ASTAR,BSTAR,CSTAR,ALSTAR,BESTAR,GASTAR
      L(1,1)=A*SINB*SINGST
      L(1,2)=0.
      L(1,3)=0.
      L(2,1)=-1.*A*SINB*COSGST
      L(2,2)=B*SINA
      L(2,3)=0.
      L(3,1)=A*COSB
      L(3,2)=B*COSA
      L(3,3)=C
      WRITE(12,503) ((L(I,J),J=1,3),I=1,3)
C     TRANSPOSE OF L MATRIX
      CALL MATRANS(L,LT)
C     INVERT L^T
      CALL MATINV(LT,U,D)
      WRITE(12,504) ((U(I,J),J=1,3),I=1,3)
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
C     GET INVERSE OF METRIC TENSOR
      CALL MATINV(G,GINV,D)
      WRITE(*,511)
125    WRITE(*,515) 
       READ(*,'(A4)') TEST
       IF(TEST.EQ.'GC'.OR.TEST.EQ.'gc') GO TO 126
       IF(TEST.EQ.'tp'.OR.TEST.EQ.'TP') GO TO 190
       IF(TEST.EQ.'END'.OR.TEST.EQ.'end') GO TO 250
       go to 125
126    WRITE(*,520)
       READ(*,'(A4)') TEST
       IF(TEST.EQ.'END'.OR.TEST.EQ.'end') GO TO 250
       IF(TEST.EQ.'rd'.OR.TEST.EQ.'RD') GO TO 155
       IF(TEST.EQ.'rr'.OR.TEST.EQ.'RR') GO TO 145
       IF(TEST.EQ.'DD'.OR.TEST.EQ.'dd') GO TO 135
       IF(TEST.EQ.'pr'.OR.TEST.EQ.'PR') GO TO 175
       go to 125
135   LL=1
      LIN=LL
140   WRITE(*,542)
      READ(*,*) HD(LL),KD(LL),LD(LL)
      WRITE(*,540)
      READ(*,*) HD2(LL),KD2(LL),LD2(LL)
      C1D(LL)=L(1,1)*HD(LL)+L(1,2)*KD(LL)+L(1,3)*LD(LL)
      C2D(LL)=L(2,1)*HD(LL)+L(2,2)*KD(LL)+L(2,3)*LD(LL)
      C3D(LL)=L(3,1)*HD(LL)+L(3,2)*KD(LL)+L(3,3)*LD(LL)
      C1D2(LL)=L(1,1)*HD2(LL)+L(1,2)*KD2(LL)+L(1,3)*LD2(LL)
      C2D2(LL)=L(2,1)*HD2(LL)+L(2,2)*KD2(LL)+L(2,3)*LD2(LL)
      C3D2(LL)=L(3,1)*HD2(LL)+L(3,2)*KD2(LL)+L(3,3)*LD2(LL)
      CALL DOTP(C1D(LL),C2D(LL),C3D(LL),C1D2(LL),C2D2(LL),C3D2(LL),
     1ANG(LL),DUM)
      LL=LL+1
        write(*,550)
	READ(*,'(A4)') TEST
	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
     1TEST.EQ.'no')GO TO 141
      go to 140
141   LOUT=LL-1
      WRITE(12,527)
      WRITE(*,527)
      WRITE(12,528) 
      WRITE(*,528)
      DO 142 LL=LIN,LOUT
      WRITE(12,510) HD(LL),KD(LL),LD(LL),HD2(LL),KD2(LL),LD2(LL),
     1ANG(LL)
142   WRITE(*,510) HD(LL),KD(LL),LD(LL),HD2(LL),KD2(LL),LD2(LL),ANG(LL)
      GO TO 125
145   LL=1
      LIN=LL
150	WRITE(*,530)
      READ(*,*) HR(LL),KR(LL),LR(LL)
      WRITE(*,531)
      READ(*,*) HR2(LL),KR2(LL),LR2(LL)
      C1R(LL)=U(1,1)*HR(LL)+U(1,2)*KR(LL)+U(1,3)*LR(LL)
      C2R(LL)=U(2,1)*HR(LL)+U(2,2)*KR(LL)+U(2,3)*LR(LL)
      C3R(LL)=U(3,1)*HR(LL)+U(3,2)*KR(LL)+U(3,3)*LR(LL)
      C1R2(LL)=U(1,1)*HR2(LL)+U(1,2)*KR2(LL)+U(1,3)*LR2(LL)
      C2R2(LL)=U(2,1)*HR2(LL)+U(2,2)*KR2(LL)+U(2,3)*LR2(LL)
      C3R2(LL)=U(3,1)*HR2(LL)+U(3,2)*KR2(LL)+U(3,3)*LR2(LL)
      CALL DOTP(C1R(LL),C2R(LL),C3R(LL),C1R2(LL),C2R2(LL),C3R2(LL),
     1ANG(LL),DUM)      
      LL=LL+1
        write(*,550)
	READ(*,'(A4)') TEST
	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
     1TEST.EQ.'no')GO TO 151
      go to 150
151   LOUT=LL-1
      WRITE(12,521)
      WRITE(*,521)
      WRITE(12,522) 
      WRITE(*,522)
      DO 152 LL=LIN,LOUT
      WRITE(12,510) HR(LL),KR(LL),LR(LL),HR2(LL),KR2(LL),LR2(LL),
     1ANG(LL)
152   WRITE(*,510) HR(LL),KR(LL),LR(LL),HR2(LL),KR2(LL),LR2(LL),ANG(LL)
      GO TO 125
C     OBLIQUITY CALCULATION (OR, OTHERWISE STATED,
C     ANGLE BETWEEN RECIP. VECTOR AND DIRECT SPACE PLANE)
155   LL=1
      LIN=LL
160	WRITE(*,530)
      READ(*,*) HR(LL),KR(LL),LR(LL)
      WRITE(*,540)
      READ(*,*) HD(LL),KD(LL),LD(LL)
      C1R(LL)=U(1,1)*HR(LL)+U(1,2)*KR(LL)+U(1,3)*LR(LL)
      C2R(LL)=U(2,1)*HR(LL)+U(2,2)*KR(LL)+U(2,3)*LR(LL)
      C3R(LL)=U(3,1)*HR(LL)+U(3,2)*KR(LL)+U(3,3)*LR(LL)
      C1D(LL)=L(1,1)*HD(LL)+L(1,2)*KD(LL)+L(1,3)*LD(LL)
      C2D(LL)=L(2,1)*HD(LL)+L(2,2)*KD(LL)+L(2,3)*LD(LL)
      C3D(LL)=L(3,1)*HD(LL)+L(3,2)*KD(LL)+L(3,3)*LD(LL)
      CALL DOTP(C1R(LL),C2R(LL),C3R(LL),C1D(LL),C2D(LL),C3D(LL),
     1ANG(LL),DUM)
        LL=LL+1
        write(*,550)
	READ(*,'(A4)') TEST
	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
     1TEST.EQ.'no')GO TO 165
      GO TO 160
165   LOUT=LL-1
      WRITE(12,523)
      WRITE(*,523)
      WRITE(12,524) 
      WRITE(*,524)
      DO 170 J=LIN,LOUT
      WRITE(12,510) HR(J),KR(J),LR(J),HD(J),KD(J),LD(J),ANG(J)
170   WRITE(*,510) HR(J),KR(J),LR(J),HD(J),KD(J),LD(J),ANG(J)
      GO TO 125
C
C     Calculate vector UU : the projection of vector VV on a plane 
C     with normal NN, i.e., UU = vv - (VV dot NN)NN
C
175   LL=1
      LIN=LL
180	WRITE(*,535)
      READ(*,*) HR(LL),KR(LL),LR(LL)
      WRITE(*,536)
      READ(*,*) HR2(LL),KR2(LL),LR2(LL)
      WRITE(*,541)
      READ(*,*) HD(LL),KD(LL),LD(LL)
      C1R(LL)=U(1,1)*HR(LL)+U(1,2)*KR(LL)+U(1,3)*LR(LL)
      C2R(LL)=U(2,1)*HR(LL)+U(2,2)*KR(LL)+U(2,3)*LR(LL)
      C3R(LL)=U(3,1)*HR(LL)+U(3,2)*KR(LL)+U(3,3)*LR(LL)
      C1R2(LL)=U(1,1)*HR2(LL)+U(1,2)*KR2(LL)+U(1,3)*LR2(LL)
      C2R2(LL)=U(2,1)*HR2(LL)+U(2,2)*KR2(LL)+U(2,3)*LR2(LL)
      C3R2(LL)=U(3,1)*HR2(LL)+U(3,2)*KR2(LL)+U(3,3)*LR2(LL)
      CRL=SQRT(C1R(LL)**2+C2R(LL)**2+C3R(LL)**2)
      CR2L=SQRT(C1R2(LL)**2+C2R2(LL)**2+C3R2(LL)**2)
      VV1(LL)=C1R(LL)/CRL
      VV2(LL)=C2R(LL)/CRL
      VV3(LL)=C3R(LL)/CRL
      NN1(LL)=C1R2(LL)/CR2L
      NN2(LL)=C2R2(LL)/CR2L
      NN3(LL)=C3R2(LL)/CR2L
      CALL DOTP(NN1(LL),NN2(LL),NN3(LL),VV1(LL),VV2(LL),VV3(LL),
     1ANG(LL),LEN(LL))
      UU1(LL)=VV1(LL)-LEN(LL)*NN1(LL)
      UU2(LL)=VV2(LL)-LEN(LL)*NN2(LL)
      UU3(LL)=VV3(LL)-LEN(LL)*NN3(LL)
      C1D(LL)=L(1,1)*HD(LL)+L(1,2)*KD(LL)+L(1,3)*LD(LL)
      C2D(LL)=L(2,1)*HD(LL)+L(2,2)*KD(LL)+L(2,3)*LD(LL)
      C3D(LL)=L(3,1)*HD(LL)+L(3,2)*KD(LL)+L(3,3)*LD(LL)
      CALL DOTP(UU1(LL),UU2(LL),UU3(LL),C1D(LL),C2D(LL),C3D(LL),
     1ANG(LL),LEN(LL))
      LL=LL+1
        write(*,550)
	READ(*,'(A4)') TEST
	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
     1TEST.EQ.'no')GO TO 185
      GO TO 180
185   LOUT=LL-1
      write(12,525)
      WRITE(*,525)
      WRITE(12,526)
      WRITE(*,526)
      DO 186 LL=LIN,LOUT
      WRITE(12,560) HR(LL),KR(LL),LR(LL),HR2(LL),KR2(LL),LR2(LL),
     1HD(LL),KD(LL),LD(LL),ANG(LL)
186   WRITE(*,560) HR(LL),KR(LL),LR(LL),HR2(LL),KR2(LL),LR2(LL),
     1HD(LL),KD(LL),LD(LL),ANG(LL)
      GO TO 125
190   WRITE(*,565) 
      READ(*,'(A4)') TEST
      IF(TEST.EQ.'mat'.OR.TEST.EQ.'MAT') GO TO 200
      IF(TEST.EQ.'PL'.OR.TEST.EQ.'pl') GO TO 195
      IF(TEST.EQ.'DI'.OR.TEST.EQ.'di') GO TO 191
      IF(TEST.EQ.'END'.OR.TEST.EQ.'end') go to 125
      go to 190
C     GET APPROX.[HKL] to HKL TWIN DIRECTION
C     OBTAINED USING g* x HKL (recip)
C     N. B. HKL IS COVARIANT 
191   WRITE(*,566)
      READ(*,*) HKL
      CALL MATVEC(HKL,GINV,UVWP)
      VAL1=MAXVAL(UVWP)
      DO 192 I=1,3
192   UVWPM(I)=UVWP(I)/VAL1
      WRITE(*,567) UVWP,UVWPM
      GO TO 190      
C     GET APPROX. HKL to [HKL] TWIN DIRECTION
C     OBTAINED USING g x HKL (DIRECT)
C     N.B. [HKL] IS CONTRAVARIANT
195   WRITE(*,568)
      READ(*,*) UVW
      CALL MATVEC(UVW,G,HKLP)
      VAL1=MAXVAL(HKLP)
      DO 196 I=1,3
196   HKLPM(I)=HKLP(I)/VAL1
      WRITE(*,569) HKLP,HKLPM
      GO TO 190
200   WRITE(*,532)
      READ(*,*) HR(1),KR(1),LR(1)
      IF(HR(1).EQ.100) GO TO 220
      C1R(1)=U(1,1)*HR(1)+U(1,2)*KR(1)+U(1,3)*LR(1)
      C2R(1)=U(2,1)*HR(1)+U(2,2)*KR(1)+U(2,3)*LR(1)
      C3R(1)=U(3,1)*HR(1)+U(3,2)*KR(1)+U(3,3)*LR(1)
      CRL=SQRT(C1R(1)**2+C2R(1)**2+C3R(1)**2)
      VV1(1)=C1R(1)/CRL
      VV2(1)=C2R(1)/CRL
      VV3(1)=C3R(1)/CRL
      RP(1,1)=2.*VV1(1)**2-1.
      RP(2,2)=2.*VV2(1)**2-1.
      RP(3,3)=2.*vv3(1)**2-1.
      RP(1,2)=2.*VV1(1)*VV2(1)
      RP(1,3)=2.*VV1(1)*VV3(1)
      RP(2,3)=2.*VV2(1)*VV3(1)
      RP(2,1)=RP(1,2)
      RP(3,1)=RP(1,3)
      RP(3,2)=RP(2,3)
      CALL MATMULT(RP,U,3,TEMP)
      CALL MATMULT(LT,TEMP,3,R)
      WRITE(*,570) HR(1),KR(1),LR(1)
      WRITE(12,570) HR(1),KR(1),LR(1)
      WRITE(*,571) ((R(I,J), J=1,3),I=1,3)
      WRITE(12,571) ((R(I,J), J=1,3),I=1,3)
220   WRITE(*,543)
      READ(*,*) HD(1),KD(1),LD(1)
      IF(HD(1).EQ.100) GO TO 190
      C1D(1)=L(1,1)*HD(1)+L(1,2)*KD(1)+L(1,3)*LD(1)
      C2D(1)=L(2,1)*HD(1)+L(2,2)*KD(1)+L(2,3)*LD(1)
      C3D(1)=L(3,1)*HD(1)+L(3,2)*KD(1)+L(3,3)*LD(1)
      CRD=SQRT(C1D(1)**2+C2D(1)**2+C3D(1)**2)
      VV1D(1)=C1D(1)/CRD
      VV2D(1)=C2D(1)/CRD
      VV3D(1)=C3D(1)/CRD
      RPD(1,1)=2.*VV1D(1)**2-1.
      RPD(2,2)=2.*VV2D(1)**2-1.
      RPD(3,3)=2.*vv3D(1)**2-1.
      RPD(1,2)=2.*VV1D(1)*VV2D(1)
      RPD(1,3)=2.*VV1D(1)*VV3D(1)
      RPD(2,3)=2.*VV2D(1)*VV3D(1)
      RPD(2,1)=RPD(1,2)
      RPD(3,1)=RPD(1,3)
      RPD(3,2)=RPD(2,3)
      CALL MATMULT(RPD,U,3,TEMP)
      CALL MATMULT(LT,TEMP,3,RD)
      WRITE(*,572) HD(1),KD(1),LD(1)
      WRITE(12,572) HD(1),KD(1),LD(1)
      WRITE(*,571) ((RD(I,J), J=1,3),I=1,3)
      WRITE(12,571) ((RD(I,J), J=1,3),I=1,3)
      GO TO 190
C      
500   FORMAT(//' CELLTOOL 2.0   BY BRUCE FOXMAN    MAY 2021'//
     1' ENTER A, B, C, ALPHA, BETA, GAMMA : ',$)
501   FORMAT(//' DIRECT CELL CONSTANTS: ', 6F10.4)
502   FORMAT(//' RECIPROCAL CELL CONSTANTS:' ,3F10.6, 3F10.4)
503   FORMAT(//' DIRECT ORTHOGONALISATION MATRIX, J. S. ROLLETT,'
     1' COMPUTING METHODS IN CRYSTALLOGRAPHY,'/' PERGAMON, 1965, PP. '
     2,'22-24: ',//3(3F12.4,/))
504   FORMAT(//' RECIP. ORTHOGONALISATION MATRIX, J. S. ROLLETT,'
     1' COMPUTING METHODS IN CRYSTALLOGRAPHY,'/' PERGAMON, 1965, PP. '     
     2,'22-24: ',//3(3F12.6,/))
510   FORMAT(/6(F7.2,2x),F10.3)
511   FORMAT(//' Note that, for the options listed below: '//,
     1' h1 k1 l1 ... h2 k2 l2 etc. [no brackets] refer to a r.l. ',
     2'vector'
     3//' [h1 k1 l1] etc. , with brackets, refer to a direct space '
     4'direction',//' (h1 k1 l1) etc., with parentheses, refer to '
     5'planes'//)
515   FORMAT(/' Enter [TP] to calculate TWIN PARAMETERS',/
     1'       [GC] to calculate OR GENERAL CELL ANGLES',/
     2' OR    [END] TO STOP PROGRAM '//' Choice: - ',$)
520   FORMAT(/' Options CALCULATE : '
     1,//' [RD] get h1 k1 l1 to [h2 k2 l2] ANGLE (reciprocal '
     2'to direct angle (= obliquity))'
     3 //, ' [RD] get h1 k1 l1 to h2 k2 l2 ANGLE (reciprocal '
     4,'to reciprocal angle) '
     5 //, ' [DD] get [h1 k1 l1] to [h2 k2 l2] ANGLE (direct to'
     6,' direct angle)'
     7 //, ' [PR] get angle between h1 k1 l1 projection onto plane',
     8' (h2 k2 l2) and [h3 k3 l3]',
     9 //' ( OR END )',//' ENTER [RD, RR, DD, PR, END]: ',$)
521   FORMAT(//' CALCULATION OF ANGLE BETWEEN TWO R.L. VECTORS ',
     1'OR TWO DIRECT SPACE PLANES')
522   FORMAT(//'    h1       k1       l1     ',
     1'  h2       k2       l2        ANGLE')
523   FORMAT(//' CALCULATION OF ANGLE BETWEEN A R.L. VECTOR ',
     1'OR DIRECT SPACE PLANE'/' AND A DIRECT SPACE DIRECTION',//,
     2'  THIS IS ALSO EQUIVALENT TO AN OBLIQUITY CALCULATION')
524   FORMAT(//'    h1       k1       l1    ',
     1'  [h2       k2       l2]      ANGLE      (obliquity)')
525   FORMAT(//' CALCULATION OF THE PROJECTION',
     1' OF A R.L. VECTOR h1 k1 l1 ONTO ',/' A DIRECT SPACE PLANE ',
     2'(h2 k2 l2) AND THE ANGLE BETWEEN THE ',/' PROJECTION AND A '
     3,'DIRECT LATTICE DIRECTION [h3 k3 l3]')
526   FORMAT(//'    h1         k1       l1        (h2       k2        '
     1'l2)       [h3       k3       l3]        ANGLE')
527   FORMAT(//' CALCULATION OF ANGLE BETWEEN TWO DIRECT SPACE',
     1' DIRECTIONS')
528   FORMAT(//' [h1     k1     l1]   [h2     k2     l2]        ANGLE')
530   FORMAT (//' ENTER h1 k1 l1 r.l. vector (or (h1 k1 l1)): ',$)
531   FORMAT (//' ENTER h2 k2 l2 r.l. vector (or (h2 k2 l2)): ',$)
532   FORMAT (//' ENTER h1 k1 l1 r.l. vector (h = 100 to skip to ',
     1'[h k l] calc.: ',$)
535   FORMAT (//' ENTER h1 k1 l1 r.l. vector: ',$)
536   FORMAT (//' ENTER (h2 k2 l2) plane normal: ',$)
540   FORMAT (//' ENTER [h2 k2 l2] direction: ',$)
541   FORMAT (//' ENTER [h3 k3 l3] direction: ',$)
542   FORMAT (//' ENTER [h1 k1 l1] direction: ',$)
543   FORMAT (//' ENTER [h1 k1 l1] direction (h = 100 to skip to menu:'
     1,' ',$)
550   FORMAT (//' ANOTHER? [Y OR N]: ',$)
560   FORMAT(/6(F7.2,3X),3(F7.2,3X),F10.3)
565   FORMAT(//' PREDICTION OF DIRECTION NORMAL TO LATTICE PLANE',/,
     1' OR INDICES OF PLANE NORMAL TO LATTICE ROW.',//,
     2' ENTER [DI] TO PREDICT DIRECTION, [PL] TO PREDICT PLANE OR,',/,
     3' [MAT] TO DERIVE TWIN LAW MATRIX FROM SPECIFIC 180 DEG ROTATION',
     4 / ' OR [END] TO RETURN TO MAIN MENU: ',$)
566   FORMAT (//' ENTER (h k l) plane: ',$)
567   FORMAT (//' RAW [u''v''w''] VALUES :',3F10.5,//,
     1' NORMALIZED [u''v''w''] VALUES: ',3F8.3,//,
     2' VALUES ARE NON-INTEGRAL; STUDY TO CHOOSE BEST DIRECTION') 
568   FORMAT (//' ENTER [u v w] direction: ',$)
569   FORMAT (//' RAW (h''k''l'') VALUES :',3F8.3,//,
     1' NORMALIZED (h''k''l'') VALUES: ',3F8.3,//
     2' VALUES ARE NON-INTEGRAL; STUDY TO CHOOSE BEST PLANE') 
570   FORMAT(//' TWIN LAW FOR 180 DEG ROTATION ABOUT R.L. VECTOR: ',
     13F8.3)
571   FORMAT(//,3(3F10.6,/))
572   FORMAT(//' TWIN LAW FOR 180 DEG ROTATION ABOUT DIRECT AXIS: [',
     13F8.3, ']')
250   STOP
      END
      SUBROUTINE MATRANS(A, AT)
      implicit real*8 (a-h,o-z)
      REAL*8 A(3,3), AT(3,3)
      DO 10 I = 1, 3
      DO 20 J = 1, 3
         AT(J,I) = A(I,J)
20    CONTINUE
10    CONTINUE
      RETURN
      END
      SUBROUTINE MATINV(A,B,D)
      implicit real*8 (a-h,o-z)
      INTEGER*2 K,N
c      REAL*4    A,B,D
      DIMENSION A(3,3),B(3,3)
      CALL VMULT(A(1,2),A(1,3),B(1,1))
      CALL VMULT(A(1,3),A(1,1),B(1,2))
      CALL VMULT(A(1,1),A(1,2),B(1,3))
      D=A(1,1)*B(1,1)+A(2,1)*B(2,1)+A(3,1)*B(3,1)
      IF(D.EQ.0.0)RETURN
      N=1
10    K=1
15    B(K,N)=B(K,N)/D
      K=K+1
      IF(K.LE.3)GOTO 15
      N=N+1
      IF(N.LE.3)GOTO 10
      CALL MACOL(B)
      RETURN
      END
      SUBROUTINE DOTP(A1,A2,A3,B1,B2,B3,ANG,CANG)
      implicit real*8 (a-h,o-z)
      ANORM=SQRT(A1**2+A2**2+A3**2)
      BNORM=SQRT(B1**2+B2**2+B3**2)
      AN1=A1/ANORM
      AN2=A2/ANORM
      AN3=A3/ANORM
      BN1=B1/BNORM
      BN2=B2/BNORM
      BN3=B3/BNORM
      CANG=AN1*BN1+AN2*BN2+AN3*BN3
      ANG=57.295*ACOS(CANG)
      RETURN
      END
      subroutine matmult(a,b,n,c)
      REAL*8 a(n,n),b(n,n),c(n,n)
      do 100 i=1,n
      do 100 j=1,n
      c(i,j)=0.
      do 100 k=1,n
  100 c(i,j)=c(i,j)+a(i,k)*b(k,j)
      return
      end
      SUBROUTINE VMULT(A,B,C)
      implicit real*8 (a-h,o-z)
      DIMENSION A(3),B(3),C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
      SUBROUTINE MACOL(A)
      implicit real*8 (a-h,o-z)
      INTEGER*2 K,N
C     REAL*4    A,T
      DIMENSION A(3,3)
      N=1
 5    K=2
10    T=A(K,N)
      A(K,N)=A(N,K)
      A(N,K)=T
      IF(N.EQ.3)RETURN
      K=K+1
      IF(K.LE.3)GOTO 10
      N=3
      GOTO 5
      END
      SUBROUTINE MATVEC(A,B,V)
C
      implicit real*8 (a-h,o-z)
      DIMENSION A(3,1),B(3,3),C(3,1),V(3)
C       A = VECTOR, B = 3 X 3, C = RESULT
C
      DO 10 I=1,3
      V(I)=0. 
      DO 10 J=1,3
10    V(I)=V(I)+B(J,I)*A(J,1) 
C       VMOD=SQRT(V(1)**2+V(2)**2+V(3)**2)  
C       DO 20 I=1,3
C 20      C(I,1)=V(I)/VMOD
      RETURN
      END
