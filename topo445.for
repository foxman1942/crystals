C     PROGRAM TOPO16a 4.45 04.03.20
      DIMENSION RM(3,3),RD(3,3),RINVM(3,3),RINVD(3,3),TM(3,3),
     1THETAM(50),CHIM(50),PHIM(50),THETAD(50),CHID(5),ANG(9),
     2PHID(50),C1M(50),C2M(50),C3M(50),C1D(50),C2D(50),C3D(50),
     3TMINV(3,3),CELL(7),CELLM(7),CELLD(7)
      DIMENSION RINVTM(3,3),RINVTD(3,3),RANG(9),RC1M(50),RC2M(50),
     1RC3M(50),RC1D(50),RC2D(50),RC3D(50),RINVTMI(3,3),RINVTDI(3,3),
     2R1M(50),R2M(50),R3M(50),D1M(50),D2M(50),D3M(50),
     3R1D(50),R2D(50),R3D(50),D1D(50),D2D(50),D3D(50)
      REAL HM(50),KM(50),LM(50),HD(50),KD(50),LD(50),
     1HMR(50),HMD(50),KMR(50),KMD(50),LMR(50),LMD(50),
     2HDR(50),HDD(50),KDR(50),KDD(50),LDR(50),LDD(50)
      COMMON /DATC/CELL
      CHARACTER*32 FILMOTH, FILDAU, FILOUT, CR, DATTYPE
      CHARACTER*4 TEST
      DATA HM(1),HM(2),HM(3),KM(1),KM(2),KM(3),LM(1),LM(2),LM(3)/
     11.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/
      DATA HD(1),HD(2),HD(3),KD(1),KD(2),KD(3),LD(1),LD(2),LD(3)/
     11.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/
      DATA NDATA,NEXPERT,NTWIN,NMOTH,NDAU/0,0,0,1,1/
C
90    WRITE(*,300)
      READ(*,'(A4)') TEST
      IF(TEST.EQ.'Y'.OR.TEST.EQ.'y') NEXPERT=1
      IF(TEST.EQ.'') GO TO 90
      WRITE(*,301)
      READ(*,'(A32)') DATTYPE
      IF(DATTYPE.EQ.'CIF'.OR.DATTYPE.EQ.'cif') NDATA=1
C     IF(DATTYPE.EQ.'') 
      IF(NDATA.EQ.1) GO TO 92
91    WRITE(*,302)
      READ(*,'(A4)') TEST
      IF(TEST.EQ.'Y'.OR.TEST.EQ.'y') NTWIN=1
      IF(TEST.EQ.'') GO TO 91
92    WRITE(*,305)
      READ(*,'(A32)') FILMOTH
      IF (NTWIN.EQ.0) GO TO 93
      WRITE(*,306)
      READ(*,*) NMOTH
93    WRITE(*,310)
      READ(*,'(A32)') FILDAU
      IF (NTWIN.EQ.0) GO TO 94
      WRITE(*,311)
      READ(*,*) NDAU
94    WRITE(*,312)
      READ(*,'(A32)') FILOUT
      IF(FILOUT.EQ.'') FILOUT='TOPO.OUT'
      OPEN(10,STATUS='OLD',FILE=FILMOTH,ACCESS='SEQUENTIAL',
     1FORM='FORMATTED')
      OPEN(11,STATUS='OLD',FILE=FILDAU,ACCESS='SEQUENTIAL',
     1FORM='FORMATTED')
      OPEN(12,STATUS='UNKNOWN',FILE=FILOUT,ACCESS='SEQUENTIAL',
     1FORM='FORMATTED')
      IF(NDATA.EQ.1) GO TO 100
      READ (10,320)
      IF (NMOTH.EQ.1) GO TO 96
      DO 95 I=2,NMOTH
95    READ(10,321)
96    CALL PARSEP4P(10,CELLM)
C     READ (10,325)  CELLM
      WRITE(*,326) CELLM
      READ (10,330) ((RM(I,J),J=1,3),I=1,3)
      READ (11,320)
      write(12,326) cellm
      IF(NDAU.EQ.1) GO TO 98
      DO 97 I=2,NDAU
97    READ(11,321)
98    CALL PARSEP4P(11,CELLD)
C     READ (11,325)  CELLD
      READ (11,330) ((RD(I,J),J=1,3),I=1,3)
      WRITE(12,340) FILMOTH
      WRITE(12,341) ((RM(I,J),J=1,3),I=1,3)
C     WRITE(*,340) FILMOTH
C     WRITE(*,341) ((RM(I,J),J=1,3),I=1,3)
      CALL MATINV(RM,RINVM,D)
      WRITE(12,342) 
      WRITE(12,341) ((RINVM(I,J),J=1,3),I=1,3)
      WRITE(12,327) CELLD
      WRITE(12,343) FILDAU
      WRITE(12,341) ((RD(I,J),J=1,3),I=1,3)
C     WRITE(*,342) 
C     WRITE(*,341) ((RINVM(I,J),J=1,3),I=1,3)
C     WRITE(*,343) FILDAU
      write(*,327) CELLD
C     WRITE(*,341) ((RD(I,J),J=1,3),I=1,3)
      CALL MATINV(RD,RINVD,D)
      WRITE(12,344) 
      WRITE(12,341) ((RINVD(I,J),J=1,3),I=1,3)
C     WRITE(*,344) 
C     WRITE(*,341) ((RINVD(I,J),J=1,3),I=1,3) 
      GO TO 105
100   CALL PARSECIF(FILMOTH,10,RM)
      CALL PARSECIF(FILDAU,11,RD)
      WRITE(12,340) FILMOTH
      WRITE(12,341) ((RM(I,J),J=1,3),I=1,3)
C     WRITE(*,340) FILMOTH
C     WRITE(*,341) ((RM(I,J),J=1,3),I=1,3)
      CALL MATINV(RM,RINVM,D)
      WRITE(12,342) 
      WRITE(12,341) ((RINVM(I,J),J=1,3),I=1,3)
      WRITE(12,343) FILDAU
      WRITE(12,341) ((RD(I,J),J=1,3),I=1,3)
C     WRITE(*,342) 
C     WRITE(*,341) ((RINVM(I,J),J=1,3),I=1,3)
C     WRITE(*,343) FILDAU
C     WRITE(*,341) ((RD(I,J),J=1,3),I=1,3)
      CALL MATINV(RD,RINVD,D)
      WRITE(12,344) 
      WRITE(12,341) ((RINVD(I,J),J=1,3),I=1,3)
C     WRITE(*,344) 
C     WRITE(*,341) ((RINVD(I,J),J=1,3),I=1,3) 
105   CALL MATMULT(RINVD, RM, 3, TM)
      WRITE(*,345)
      WRITE(*,341) ((TM(I,J),J=1,3),I=1,3)
      WRITE(12,345)
      WRITE(12,341) ((TM(I,J),J=1,3),I=1,3)
      CALL MATINV(TM, TMINV, D)
      WRITE(*,346) D
      WRITE(12,346) D
      WRITE(*,347)
      WRITE(*,341) ((TMINV(I,J),J=1,3),I=1,3)
      WRITE(12,347)
      WRITE(12,341) ((TMINV(I,J),J=1,3),I=1,3)
      DO 110 JJ=1,3
      C1M(JJ)=RM(1,1)*HM(JJ)+RM(1,2)*KM(JJ)+RM(1,3)*LM(JJ)
      C2M(JJ)=RM(2,1)*HM(JJ)+RM(2,2)*KM(JJ)+RM(2,3)*LM(JJ)
      C3M(JJ)=RM(3,1)*HM(JJ)+RM(3,2)*KM(JJ)+RM(3,3)*LM(JJ)
C     WRITE(12,350) HM(JJ),KM(JJ),LM(JJ),C1M(JJ),C2M(JJ),C3M(JJ)
110   CONTINUE
      DO 120 JJ=1,3
      C1D(JJ)=RD(1,1)*HD(JJ)+RD(1,2)*KD(JJ)+RD(1,3)*LD(JJ)
      C2D(JJ)=RD(2,1)*HD(JJ)+RD(2,2)*KD(JJ)+RD(2,3)*LD(JJ)
      C3D(JJ)=RD(3,1)*HD(JJ)+RD(3,2)*KD(JJ)+RD(3,3)*LD(JJ)
C     WRITE(12,350) HD(JJ),KD(JJ),LD(JJ),C1D(JJ),C2D(JJ),C3D(JJ)
120   CONTINUE
      L=1
      WRITE(12,351)
      WRITE(*,351)
      DO 130 I=1,3
      DO 130 J=1,3
      CALL DOTP(C1M(I),C2M(I),C3M(I),C1D(J),C2D(J),C3D(J),ANG(L))
      WRITE(12,360) HM(I),KM(I),LM(I),HD(J),KD(J),LD(J),ANG(L)
      WRITE(*,360) HM(I),KM(I),LM(I),HD(J),KD(J),LD(J),ANG(L)
130   L=L+1
      WRITE(*,353)
      READ(*,'(A1)') CR
C     Test here for 180's; reset matrices if appropriate and return/stay here after recalc
C     IF(ANG(1).GE.170..AND.ANG(5).GE.170.) ITEST=1
C     IF(ANG(1).GE.170..AND.ANG(9).GE.170.) ITEST=2     
C     IF(ANG(5).GE.170..AND.ANG(9).GE.170.) ITEST=3
C     IF(ITEST.EQ.0) GO TO 134
      IF(NEXPERT.EQ.0) GO TO 134
      WRITE(*,361) 
      WRITE(*,362)
      READ(*,'(I10)') ITEST
      IF(ITEST.EQ.0) GO TO 134 
      CALL ROT(RM,ITEST)
      IF(ITEST.GT.O) CALL MATINV(RM,RINVM,D) 
C     copycat redo section
      WRITE(12,348)
      WRITE(12,363)
      WRITE(12,340) FILMOTH
      WRITE(12,341) ((RM(I,J),J=1,3),I=1,3)
      WRITE(*,363)
      WRITE(*,340) FILMOTH
      WRITE(*,341) ((RM(I,J),J=1,3),I=1,3)
      CALL MATINV(RM,RINVM,D)
      WRITE(12,363)
      WRITE(12,342) 
      WRITE(12,341) ((RINVM(I,J),J=1,3),I=1,3)
      WRITE(12,343) FILDAU
      WRITE(12,341) ((RD(I,J),J=1,3),I=1,3)
      WRITE(*,363)
      WRITE(*,342) 
      WRITE(*,341) ((RINVM(I,J),J=1,3),I=1,3)
C     WRITE(*,343) FILDAU
C     WRITE(*,341) ((RD(I,J),J=1,3),I=1,3)
      WRITE(12,344) 
      WRITE(12,341) ((RINVD(I,J),J=1,3),I=1,3)
C     WRITE(*,344) 
C     WRITE(*,341) ((RINVD(I,J),J=1,3),I=1,3)
      CALL MATMULT(RINVD, RM, 3, TM)
      WRITE(*,363)
      WRITE(*,345)
      WRITE(*,341) ((TM(I,J),J=1,3),I=1,3)
      WRITE(12,363)
      WRITE(12,345)
      WRITE(12,341) ((TM(I,J),J=1,3),I=1,3)
      CALL MATINV(TM, TMINV, D)
      WRITE(*,363)
      WRITE(*,347)
      WRITE(*,341) ((TMINV(I,J),J=1,3),I=1,3)
      WRITE(12,363)
      WRITE(12,347)
      WRITE(12,341) ((TMINV(I,J),J=1,3),I=1,3)
      DO 1101 JJ=1,3
      C1M(JJ)=RM(1,1)*HM(JJ)+RM(1,2)*KM(JJ)+RM(1,3)*LM(JJ)
      C2M(JJ)=RM(2,1)*HM(JJ)+RM(2,2)*KM(JJ)+RM(2,3)*LM(JJ)
      C3M(JJ)=RM(3,1)*HM(JJ)+RM(3,2)*KM(JJ)+RM(3,3)*LM(JJ)
C     WRITE(12,350) HM(JJ),KM(JJ),LM(JJ),C1M(JJ),C2M(JJ),C3M(JJ)
1101  CONTINUE
      DO 1201 JJ=1,3
      C1D(JJ)=RD(1,1)*HD(JJ)+RD(1,2)*KD(JJ)+RD(1,3)*LD(JJ)
      C2D(JJ)=RD(2,1)*HD(JJ)+RD(2,2)*KD(JJ)+RD(2,3)*LD(JJ)
      C3D(JJ)=RD(3,1)*HD(JJ)+RD(3,2)*KD(JJ)+RD(3,3)*LD(JJ)
C     WRITE(12,350) HD(JJ),KD(JJ),LD(JJ),C1D(JJ),C2D(JJ),C3D(JJ)
1201  CONTINUE
      L=1
      WRITE(12,351)
      WRITE(*,351)
      DO 1301 I=1,3
      DO 1301 J=1,3
      CALL DOTP(C1M(I),C2M(I),C3M(I),C1D(J),C2D(J),C3D(J),ANG(L))
      WRITE(12,360) HM(I),KM(I),LM(I),HD(J),KD(J),LD(J),ANG(L)
      WRITE(*,360) HM(I),KM(I),LM(I),HD(J),KD(J),LD(J),ANG(L)
1301  L=L+1
134   L=4
135   WRITE(*,370)
      READ(*,'(A4)') TEST
      IF(TEST.EQ.'DIRE'.OR.TEST.EQ.'dire'.OR.
     1TEST.EQ.'MD'.OR.TEST.EQ.'md'.OR.
     2TEST.EQ.'DM'.OR.TEST.EQ.'dm'.OR.
     3TEST.EQ.'IV'.OR.TEST.EQ.'iv') GO TO 136
      GO TO 135
136   IF(TEST.EQ.'DIRE'.OR.TEST.EQ.'dire') GO TO 200
      IF(TEST.EQ.'IV'.OR.TEST.EQ.'iv') GO TO 155
      IF(TEST.EQ.'DM'.OR.TEST.EQ.'dm') GO TO 145
      IF(TEST.EQ.'') GO TO 136
	WRITE(12,371)
	WRITE(12,381)
      WRITE(*,371)
140   WRITE(*,380) 
      READ (*,*) HM(L),KM(L),LM(L)
      IF(HM(L).EQ.100) GO TO 135
      WRITE(*,381)
      C1M(L)=RM(1,1)*HM(L)+RM(1,2)*KM(L)+RM(1,3)*LM(L)
      C2M(L)=RM(2,1)*HM(L)+RM(2,2)*KM(L)+RM(2,3)*LM(L)
      C3M(L)=RM(3,1)*HM(L)+RM(3,2)*KM(L)+RM(3,3)*LM(L)
      HD(L)=RINVD(1,1)*C1M(L)+RINVD(1,2)*C2M(L)+RINVD(1,3)*C3M(L)
      KD(L)=RINVD(2,1)*C1M(L)+RINVD(2,2)*C2M(L)+RINVD(2,3)*C3M(L)
      LD(L)=RINVD(3,1)*C1M(L)+RINVD(3,2)*C2M(L)+RINVD(3,3)*C3M(L)
      WRITE (12,400) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L)
      WRITE (*,400) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L)
      L=L+1
      GO TO 140 
145	WRITE(12,372)
      WRITE(12,382)
   	WRITE(*,372)
150   WRITE(*,390) 
      READ(*,*) HD(L),KD(L),LD(L)
      IF(HD(L).EQ.100) GO TO 135
      WRITE(*,382)
      C1D(L)=RD(1,1)*HD(L)+RD(1,2)*KD(L)+RD(1,3)*LD(L)
      C2D(L)=RD(2,1)*HD(L)+RD(2,2)*KD(L)+RD(2,3)*LD(L)
      C3D(L)=RD(3,1)*HD(L)+RD(3,2)*KD(L)+RD(3,3)*LD(L)
      HM(L)=RINVM(1,1)*C1D(L)+RINVM(1,2)*C2D(L)+RINVM(1,3)*C3D(L)
      KM(L)=RINVM(2,1)*C1D(L)+RINVM(2,2)*C2D(L)+RINVM(2,3)*C3D(L)
      LM(L)=RINVM(3,1)*C1D(L)+RINVM(3,2)*C2D(L)+RINVM(3,3)*C3D(L)
      WRITE (12,400) HD(L),KD(L),LD(L),HM(L),KM(L),LM(L)
      WRITE(*,400) HD(L),KD(L),LD(L),HM(L),KM(L),LM(L)
      L=L+1
      GO TO 150 
155   WRITE(12,401)
      WRITE(*,401)
      WRITE(12,402) 
      L=1
160	WRITE(*,410)
      READ(*,*) HM(L),KM(L),LM(L)
      WRITE(*,420)
      READ(*,*) HD(L),KD(L),LD(L)
      C1M(L)=RM(1,1)*HM(L)+RM(1,2)*KM(L)+RM(1,3)*LM(L)
      C2M(L)=RM(2,1)*HM(L)+RM(2,2)*KM(L)+RM(2,3)*LM(L)
      C3M(L)=RM(3,1)*HM(L)+RM(3,2)*KM(L)+RM(3,3)*LM(L)
      C1D(L)=RD(1,1)*HD(L)+RD(1,2)*KD(L)+RD(1,3)*LD(L)
      C2D(L)=RD(2,1)*HD(L)+RD(2,2)*KD(L)+RD(2,3)*LD(L)
      C3D(L)=RD(3,1)*HD(L)+RD(3,2)*KD(L)+RD(3,3)*LD(L)
      CALL DOTP(C1M(L),C2M(L),C3M(L),C1D(L),C2D(L),C3D(L),ANG(L))
      WRITE(12,360) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L),ANG(L)
C	WRITE(*,430)
C	READ(*,'(A4)') TEST
C	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
C     1TEST.EQ.'no')GO TO 165
C      GO TO 160
165   LOUT=L-1
      WRITE(*,402)
C      DO 170 LL=LIN,LOUT
170   WRITE(*,360) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L),ANG(L)
      L=L+1
	WRITE(*,430)
	READ(*,'(A4)') TEST
	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
     1TEST.EQ.'no')GO TO 135
      GO TO 160
200   DO 201 II=1,3
      HM(II)=0.0
      KM(II)=0.0
      LM(II)=0.0
      HD(II)=0.0
      KD(II)=0.0
201   LD(II)=0.0
      HM(1)=1.0
      KM(2)=1.0
      LM(3)=1.0
      HD(1)=1.0
      KD(2)=1.0
      LD(3)=1.0
      CALL MATRANS(RINVM,RINVTM)
      CALL MATRANS(RINVD,RINVTD)
      CALL MATINV(RINVTM,RINVTMI,D)
      CALL MATINV(RINVTD,RINVTDI,D)
      DO 210 JJ=1,3
      RC1M(JJ)=RINVTM(1,1)*HM(JJ)+RINVTM(1,2)*KM(JJ)+RINVTM(1,3)*LM(JJ)
      RC2M(JJ)=RINVTM(2,1)*HM(JJ)+RINVTM(2,2)*KM(JJ)+RINVTM(2,3)*LM(JJ)
      RC3M(JJ)=RINVTM(3,1)*HM(JJ)+RINVTM(3,2)*KM(JJ)+RINVTM(3,3)*LM(JJ)
C     WRITE(12,350) HM(JJ),KM(JJ),LM(JJ),RC1M(JJ),RC2M(JJ),RC3M(JJ)
210   CONTINUE
      DO 220 JJ=1,3
      RC1D(JJ)=RINVTD(1,1)*HD(JJ)+RINVTD(1,2)*KD(JJ)+RINVTD(1,3)*LD(JJ)
      RC2D(JJ)=RINVTD(2,1)*HD(JJ)+RINVTD(2,2)*KD(JJ)+RINVTD(2,3)*LD(JJ)
      RC3D(JJ)=RINVTD(3,1)*HD(JJ)+RINVTD(3,2)*KD(JJ)+RINVTD(3,3)*LD(JJ)
C     WRITE(12,350) HD(JJ),KD(JJ),LD(JJ),RC1D(JJ),RC2D(JJ),RC3D(JJ)
220   CONTINUE
      L=1
      WRITE(12,352)
      WRITE(*,352)
      DO 230 I=1,3
      DO 230 J=1,3
      CALL DOTP(RC1M(I),RC2M(I),RC3M(I),RC1D(J),RC2D(J),RC3D(J),RANG(L))
      WRITE(12,360) HM(I),KM(I),LM(I),HD(J),KD(J),LD(J),RANG(L)
      WRITE(*,360) HM(I),KM(I),LM(I),HD(J),KD(J),LD(J),RANG(L)
230   L=L+1
      WRITE(*,353)
      READ(*,'(A1)') CR
      L=4
235   WRITE(*,375)
      READ(*,'(A4)') TEST
      IF(TEST.EQ.'OBLI'.OR.TEST.EQ.'obli'.OR.
     1TEST.EQ.'MD'.OR.TEST.EQ.'md'.OR.
     2TEST.EQ.'DM'.OR.TEST.EQ.'dm'.OR.
     3TEST.EQ.'IV'.OR.TEST.EQ.'iv') GO TO 236
      GO TO 235
236   IF(TEST.EQ.'OBLI'.OR.TEST.EQ.'obli') GO TO 271
      IF(TEST.EQ.'IV'.OR.TEST.EQ.'iv') GO TO 255
      IF(TEST.EQ.'DM'.OR.TEST.EQ.'dm') GO TO 245
      IF(TEST.EQ.'')GO TO 236
	WRITE(12,373)
	WRITE(12,383)
      WRITE(*,373)
240   WRITE(*,385) 
      READ (*,*) HM(L),KM(L),LM(L)
      IF(HM(L).EQ.100) GO TO 235
      WRITE(*,383)
      RC1M(L)=RINVTM(1,1)*HM(L)+RINVTM(1,2)*KM(L)+RINVTM(1,3)*LM(L)
      RC2M(L)=RINVTM(2,1)*HM(L)+RINVTM(2,2)*KM(L)+RINVTM(2,3)*LM(L)
      RC3M(L)=RINVTM(3,1)*HM(L)+RINVTM(3,2)*KM(L)+RINVTM(3,3)*LM(L)
      HD(L)=RINVTDI(1,1)*RC1M(L)+RINVTDI(1,2)*RC2M(L)+RINVTDI(1,3)*
     1RC3M(L)
      KD(L)=RINVTDI(2,1)*RC1M(L)+RINVTDI(2,2)*RC2M(L)+RINVTDI(2,3)*
     1RC3M(L)
      LD(L)=RINVTDI(3,1)*RC1M(L)+RINVTDI(3,2)*RC2M(L)+RINVTDI(3,3)*
     1RC3M(L)
      WRITE (12,400) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L)
      WRITE (*,400) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L)
      L=L+1
      GO TO 240 
245	WRITE(12,374)
      WRITE(12,382)
   	WRITE(*,374)
250   WRITE(*,395) 
      READ(*,*) HD(L),KD(L),LD(L)
      IF(HD(L).EQ.100) GO TO 235
      WRITE(*,382)
      RC1D(L)=RINVTD(1,1)*HD(L)+RINVTD(1,2)*KD(L)+RINVTD(1,3)*LD(L)
      RC2D(L)=RINVTD(2,1)*HD(L)+RINVTD(2,2)*KD(L)+RINVTD(2,3)*LD(L)
      RC3D(L)=RINVTD(3,1)*HD(L)+RINVTD(3,2)*KD(L)+RINVTD(3,3)*LD(L)
      HM(L)=RINVTMI(1,1)*RC1D(L)+RINVTMI(1,2)*RC2D(L)+RINVTMI(1,3)*
     1RC3D(L)
      KM(L)=RINVTMI(2,1)*RC1D(L)+RINVTMI(2,2)*RC2D(L)+RINVTMI(2,3)*
     1RC3D(L)
      LM(L)=RINVTMI(3,1)*RC1D(L)+RINVTMI(3,2)*RC2D(L)+RINVTMI(3,3)*
     1RC3D(L)
      WRITE (12,400) HD(L),KD(L),LD(L),HM(L),KM(L),LM(L)
      WRITE(*,400) HD(L),KD(L),LD(L),HM(L),KM(L),LM(L)
      L=L+1
      GO TO 250 
255   WRITE(12,403)
      WRITE(*,403)
      WRITE(12,404) 
C     LIN=L
      L=1
260	WRITE(*,411)
      READ(*,*) HM(L),KM(L),LM(L)
      WRITE(*,421)
      READ(*,*) HD(L),KD(L),LD(L)
      RC1M(L)=RINVTM(1,1)*HM(L)+RINVTM(1,2)*KM(L)+RINVTM(1,3)*LM(L)
      RC2M(L)=RINVTM(2,1)*HM(L)+RINVTM(2,2)*KM(L)+RINVTM(2,3)*LM(L)
      RC3M(L)=RINVTM(3,1)*HM(L)+RINVTM(3,2)*KM(L)+RINVTM(3,3)*LM(L)
      RC1D(L)=RINVTD(1,1)*HD(L)+RINVTD(1,2)*KD(L)+RINVTD(1,3)*LD(L)
      RC2D(L)=RINVTD(2,1)*HD(L)+RINVTD(2,2)*KD(L)+RINVTD(2,3)*LD(L)
      RC3D(L)=RINVTD(3,1)*HD(L)+RINVTD(3,2)*KD(L)+RINVTD(3,3)*LD(L)
      CALL DOTP(RC1M(L),RC2M(L),RC3M(L),RC1D(L),RC2D(L),RC3D(L),
     1RANG(L))
      WRITE(12,360) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L),RANG(L)
C	WRITE(*,430)
C	READ(*,'(A4)') TEST
C	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
C    1TEST.EQ.'no')GO TO 265
C      GO TO 260
265   LOUT=L-1
      WRITE(*,402)
C      DO 270 LL=LIN,LOUT
270   WRITE(*,360) HM(L),KM(L),LM(L),HD(L),KD(L),LD(L),RANG(L)
      L=L+1
	WRITE(*,430)
	READ(*,'(A4)') TEST
	IF (TEST.EQ.'N'.OR.TEST.EQ.'NO'.OR.TEST.EQ.'n'.OR.
     1TEST.EQ.'no')GO TO 235
      GO TO 260
271   WRITE(*,500)
      L=1
      READ(*,'(A4)') TEST
      IF(TEST.EQ.'quit'.OR.TEST.EQ.'QUIT') GO TO 2000
      IF(TEST.EQ.'OBLD'.OR.TEST.EQ.'obld') GO TO 280
      IF(TEST.EQ.'') GO TO 271
      WRITE(12,501)
      WRITE(12,502)
      WRITE(*,501)
      WRITE(*,503)
      READ(*,*) HMR(L),KMR(L),LMR(L)
      WRITE(*,504)
      READ(*,*) HMD(L),KMD(L),LMD(L)
      WRITE(*,502)
      R1M(L)=RM(1,1)*HMR(L)+RM(1,2)*KMR(L)+RM(1,3)*LMR(L)
      R2M(L)=RM(2,1)*HMR(L)+RM(2,2)*KMR(L)+RM(2,3)*LMR(L)
      R3M(L)=RM(3,1)*HMR(L)+RM(3,2)*KMR(L)+RM(3,3)*LMR(L)
      D1M(L)=RINVTM(1,1)*HMD(L)+RINVTM(1,2)*KMD(L)+RINVTM(1,3)*LMD(L)
      D2M(L)=RINVTM(2,1)*HMD(L)+RINVTM(2,2)*KMD(L)+RINVTM(2,3)*LMD(L)
      D3M(L)=RINVTM(3,1)*HMD(L)+RINVTM(3,2)*KMD(L)+RINVTM(3,3)*LMD(L)
      CALL DOTP(R1M(L),R2M(L),R3M(L),D1M(L),D2M(L),D3M(L),
     1RANG(L))
      WRITE(12,505) HMR(L),KMR(L),LMR(L),HMD(L),KMD(L),LMD(L),RANG(L)
      WRITE(*,505) HMR(L),KMR(L),LMR(L),HMD(L),KMD(L),LMD(L),RANG(L)
      GO TO 271
280   CONTINUE
      WRITE(12,511)
      WRITE(12,512)
      WRITE(*,511)
      WRITE(*,513)
      READ(*,*) HDR(L),KDR(L),LDR(L)
      WRITE(*,514)
      READ(*,*) HDD(L),KDD(L),LDD(L)
      WRITE(*,512)
      R1D(L)=RD(1,1)*HDR(L)+RD(1,2)*KDR(L)+RD(1,3)*LDR(L)
      R2D(L)=RD(2,1)*HDR(L)+RD(2,2)*KDR(L)+RD(2,3)*LDR(L)
      R3D(L)=RD(3,1)*HDR(L)+RD(3,2)*KDR(L)+RD(3,3)*LDR(L)
      D1D(L)=RINVTD(1,1)*HDD(L)+RINVTD(1,2)*KDD(L)+RINVTD(1,3)*LDD(L)
      D2D(L)=RINVTD(2,1)*HDD(L)+RINVTD(2,2)*KDD(L)+RINVTD(2,3)*LDD(L)
      D3D(L)=RINVTD(3,1)*HDD(L)+RINVTD(3,2)*KDD(L)+RINVTD(3,3)*LDD(L)
      CALL DOTP(R1D(L),R2D(L),R3D(L),D1D(L),D2D(L),D3D(L),
     1RANG(L))
      WRITE(12,505) HDR(L),KDR(L),LDR(L),HDD(L),KDD(L),LDD(L),RANG(L)
      WRITE(*,505) HDR(L),KDR(L),LDR(L),HDD(L),KDD(L),LDD(L),RANG(L)
      GO TO 271
300   FORMAT (//' TOPOTAXY CALCULATOR   VERSION 4.45  04 MAR 2020',
     1//,' NOTE: ALL NUMERICAL INPUT MUST BE SEPARATED BY COMMAS',//,
     2' WILL THIS BE AN EXPERT MODE RUN?',//
     3,' ** RECOMMENDED ONLY IF YOU ARE FAMILIAR  **'
     4,/' ** WITH DIFFRACTOMETER'
     5,' MATRIX OPERATIONS **',//,' [ENTER Y OR N]: ',$)
301   FORMAT (//' ENTER ORIENTATION MATRIX FILE TYPE, EITHER "P4P" OR'
     1,' "CIF" ',// ' [DEFAULT(CR) = P4P] : ',$)
302   FORMAT (//' ARE THERE TWIN ORIENTATIONS TO BE CONSIDERED FOR ',
     1'MOTHER OR DAUGHTER?',//,' [ENTER Y OR N]: ',$)
305   FORMAT (//' ENTER FILE NAME, WITH EXTENSION, FOR MOTHER '
     1,'CRYSTAL: ',$)
306   FORMAT (//' ENTER SEQUENCE NUMBER OF TWIN IN P4P FILE FOR MOTHER '
     1,'CRYSTAL',//,
     2' [1 = PARENT MATRIX...NO TWIN]: ',$)
310   FORMAT (//' ENTER FILE NAME, WITH EXTENSION, FOR ','DAUGHTER '
     1,'CRYSTAL: ',$)
311   FORMAT (//' ENTER SEQUENCE NUMBER OF TWIN IN P4P FILE ',
     1'FOR DAUGHTER CRYSTAL...',//,
     2' [1 = PARENT MATRIX...NO TWIN]: ',$)
312   FORMAT (//' ENTER FILE NAME FOR OUTPUT [DEFAULT=TOPO.OUT]: ',$)
320   FORMAT(///)
321   FORMAT(//////)
325   FORMAT(6X, 6F10.5,F11.5,/)
326   FORMAT(/,' MOTHER CELL: ',3X,3F9.5,3F8.3,F11.3)
327   FORMAT(/,' DAUGHTHER CELL: ',3F9.5,3F8.3,F11.3)
330   FORMAT(3(7X, 3E16.7,/))
340   FORMAT(//' MOTHER CRYSTAL ORIENTATION [UB]m MATRIX FROM ', A32//)
341   FORMAT(3(7X, 3F15.8,/))
342   FORMAT(//' INVERSE OF MOTHER CRYSTAL ORIENTATION MATRIX',
     1' ([UB]m)^-1'//)
343   FORMAT(//' DAUGHTER CRYSTAL ORIENTATION [UB]d MATRIX FROM '
     1, A32//)
344   FORMAT(//' INVERSE OF DAUGHTER CRYSTAL ORIENTATION MATRIX',
     1' ([UB]d)^-1'//)
345   FORMAT(//' TOPOTACTIC TRANSFORMATION MATRIX FOR DIRECT CELLS',
     1' (PHI)'/' GOUGOUTAS, J. Z. ISR. J. CHEM. 1972, 10, 395-407'
     2//)
346   FORMAT(//' DETERMINANT OF TOPOTACTIC TRANSFORMATION MATRIX:',2X,
     1 F6.3, /, ' ***THIS SHOULD EQUAL V(DAUGHTER)/V(MOTHER)***')
347   FORMAT(//' INVERSE TOPOTACTIC TRANSFORMATION MATRIX',
     1/' FOR RECIPROCAL CELLS (PHI)^-1',/
     2' GOUGOUTAS, J. Z. ISR. J. CHEM. 1972, 10, 395-407'//)
348   FORMAT(//' RESULTS AFTER EXPERT-MODE ROTATION' //)
350   FORMAT(7X,3F4.0,3F10.7)
351   FORMAT(//' Reciprocal: Part IR -- AUTOMATIC GENERATION OF ANGLES',
     1/,       '                      BETWEEN PRIMARY R. L. VECTORS',
     2//       '                      M = MOTHER, D= DAUGHTER',
     3//'  HM   KM   LM   HD   KD   LD      ANGLE')
352   FORMAT(//' Direct: PART ID -- AUTOMATIC GENERATION OF ANGLES',/,
     1         '                  BETWEEN PRIMARY DIRECT LATTICE ',
     2'VECTORS'//'                  M = MOTHER, D= DAUGHTER',
     3//' [HM   KM   LM] [HD   KD   LD]     ANGLE')
353   FORMAT(' TYPE<CR> to CONTINUE: ',$)
360   FORMAT(/6F5.0,F10.3)
361   FORMAT(/' **************  expert mode option  ************  ',//,
     1' YOU MAY ADJUST THE ORIENTATION OF THE MOTHER CRYSTAL IN ORDER',/
     2' TO OBTAIN THE BEST FIT BETWEEN MOTHER AND DAUGHTER CELLS'//)
362   FORMAT(//' ENTER 1, 2 OR 3 FOR A ROTATION OF 180 DEGREES OF THE ',
     1/' MOTHER CRYSTAL ABOUT DIFFRACTOMETER x, y, or z, RESPECTIVELY,'
     2,/,' or <CR> FOR NO ROTATION :  ',$)
363   FORMAT(//' (****** REVISED ******) ',$)
370   FORMAT(/' CALCULATE M-->D, D-->M CORRESPONDING VECTORS, ', 
     1/' INTERVECTOR ANGLES, OR GO TO DIRECT LATTICE CALC [MD,DM,IV,',
     2'DIRECT]: ',$)
371   FORMAT(//' Reciprocal: Part IIR(a) -- GENERATE DAUGHTER INDICES ',
     1'THAT CORRESPOND'/
     2         '                        TO INPUT MOTHER INDICES',//
     3         '                        NOTE THAT DAUGHTER ',
     4 'INDICES MAY BE FRACTIONAL!')
372   FORMAT(//' Reciprocal: Part IIR(b) -- GENERATE MOTHER INDICES ',
     1'THAT CORRESPOND',/,
     2         '                        TO INPUT DAUGHTER INDICES',//
     3         '                        NOTE THAT MOTHER INDICES', 
     4 ' MAY BE FRACTIONAL!')
373   FORMAT(//' Direct: Part IID(a) -- GENERATE DAUGHTER DIRECTIONS ',
     1'THAT CORRESPOND'/
     2         '                    TO INPUT MOTHER DIRECTIONS',//
     3         '                    NOTE THAT DAUGHTER ',
     4 'DIRECTIONS MAY BE FRACTIONAL!')
374   FORMAT(//' Direct: Part IID(b) -- GENERATE MOTHER DIRECTIONS THAT',
     1' CORRESPOND',/,
     2         '                    TO INPUT DAUGHTER DIRECTIONS',//
     3         '                    NOTE THAT MOTHER DIRECTIONS', 
     4 ' MAY BE FRACTIONAL!')
375   FORMAT(/' CALCULATE M-->D, D-->M CORRESPONDING VECTORS, ', 
     1/' INTERVECTOR ANGLES, OR GO TO OBLIQUITY CALC [MD,DM,IV,OBLI]: '
     2,$)
380   FORMAT (//' ENTER H, K, L FOR MOTHER CRYSTAL ',
     1'[100,0,0 TO END]: ',$)
381   FORMAT(//'              MOTHER             ',
     1         '          DAUGHTER            ',/,
     2    7X,'H',9X,'K',9X,'L',9X,'H',9X,'K',9X,'L',
     3       / '_________________________________',
     4         '______________________________')
382   FORMAT(//'             DAUGHTER           ',
     1         '            MOTHER             ',/,
     2    7X,'H',9X,'K',9X,'L',9X,'H',9X,'K',9X,'L',
     3       / '_________________________________',
     4         '______________________________')
383   FORMAT(//'              MOTHER             ',
     1         '          DAUGHTER            ',/,
     2    6X,'[H',9X,'K',9X,'L]',8X,'[H',9X,'K',8X,'L]',
     3       / '_________________________________',
     4         '______________________________')
384   FORMAT(//'             DAUGHTER           ',
     1         '            MOTHER             ',/,
     2    7X,'H',9X,'K',9X,'L',9X,'H',9X,'K',9X,'L',
     3       / '_________________________________',
     4         '______________________________')
385   FORMAT (//' ENTER [H, K, L] FOR MOTHER CRYSTAL ',
     1'[100,0,0 TO END]: ',$)
390   FORMAT (//' ENTER H, K, L FOR DAUGHTER CRYSTAL ',
     1'[100,0,0 TO END]: ',$)
395   FORMAT (//' ENTER [H, K, L] FOR DAUGHTER CRYSTAL ',
     1'[100,0,0 TO END]: ',$)
400   FORMAT(/6F10.3)
401   FORMAT(//' Reciprocal: Part IIIR -- GENERATION OF ANGLES',
     1' BETWEEN',/, 
     2         '                        SPECIFIED R. L. VECTORS',//
     3         '                        M = MOTHER, D= DAUGHTER')
402   FORMAT(//'  HM   KM   LM   HD   KD   LD      ANGLE')
403   FORMAT(//' Direct: Part IIID -- GENERATION OF ANGLES',
     1' BETWEEN',/, 
     2         '                    SPECIFIED DIRECTIONS',//
     3         '                    M = MOTHER, D= DAUGHTER')
404   FORMAT(//' [HM   KM   LM] [HD   KD   LD]     ANGLE')
410   FORMAT (//' ENTER H, K, L FOR MOTHER CRYSTAL: ',$)
411   FORMAT (//' ENTER [H, K, L] FOR MOTHER CRYSTAL: ',$)
420   FORMAT (//' ENTER H, K, L FOR DAUGHTER CRYSTAL: ',$)
421   FORMAT (//' ENTER [H, K, L] FOR DAUGHTER CRYSTAL: ',$)
430   FORMAT (//' ANOTHER? [Y OR N]: ',$)
500   FORMAT(/' CALCULATE PLANE-TO-DIRECTION ANGLE IN MOTHER [OBLM],'
     1,/,' DAUGHTER [OBLD], OR QUIT [QUIT]: ',$)
501   FORMAT(//' PART IVA- -- GENERATE MOTHER PLANE-TO-DIRECT AXIS',
     1' ANGLES ',//12X,'  (SPECIFIC OBLIQUITY CALCULATION)' )
502   FORMAT(//'   MOTHER (PLANE)      ',
     1         'MOTHER [DIRECT AXIS]     OBLIQUITY  ',//,
     2    1X,'(H',6X,'K',6X,'L)',6X,'[H',6X,'K',6X,'L]',
     3       / '_________________________________',
     4         '______________________________')
503   FORMAT (//' ENTER (H, K, L) FOR MOTHER CRYSTAL: ',$)
504   FORMAT (//' ENTER [H, K, L] FOR MOTHER CRYSTAL: ',$)
505   FORMAT(/F4.0,2F7.0,2X,3F7.0,7X,F7.3)
511   FORMAT(//' PART IVB- -- GENERATE DAUGHTER PLANE-TO-DIRECT AXIS',
     1' ANGLES ',//12X,'  (SPECIFIC OBLIQUITY CALCULATION)' )
512   FORMAT(//'  DAUGHTER (PLANE)     ',
     1         'DAUGHTER [DIRECT AXIS]     OBLIQUITY  ',//,
     2    1X,'(H',6X,'K',6X,'L)',6X,'[H',6X,'K',6X,'L]',
     3       / '_________________________________',
     4         '______________________________')
513   FORMAT (//' ENTER (H, K, L) FOR DAUGHTER CRYSTAL: ',$)
514   FORMAT (//' ENTER [H, K, L] FOR DAUGHTER CRYSTAL: ',$)
2000  STOP
      END
      SUBROUTINE MATINV(A,B,D)
      IMPLICIT  NONE
      INTEGER*2 K,N
      REAL*4    A,B,D
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
      SUBROUTINE VMULT(A,B,C)
      IMPLICIT  NONE
      REAL*4    A,B,C
      DIMENSION A(3),B(3),C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
      SUBROUTINE MACOL(A)
      IMPLICIT  NONE
      INTEGER*2 K,N
      REAL*4    A,T
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
      SUBROUTINE DOTP(A1,A2,A3,B1,B2,B3,ANG)
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
      dimension a(n,n),b(n,n),c(n,n)
      do 100 i=1,n
      do 100 j=1,n
      c(i,j)=0.
      do 100 k=1,n
  100 c(i,j)=c(i,j)+a(i,k)*b(k,j)
      return
      end
      SUBROUTINE MATRANS(A, AT)
      REAL A(3,3), AT(3,3)
      DO 10 I = 1, 3
      DO 20 J = 1, 3
         AT(J,I) = A(I,J)
20    CONTINUE
10    CONTINUE
      RETURN
      END
      subroutine rot(c,ITEST)
      dimension c(3,3)
      go to (10,20,30), ITEST
   10 c(1,2)=-c(1,2)
      c(2,2)=-c(2,2)
      c(3,2)=-c(3,2)
      c(1,3)=-c(1,3)
      c(2,3)=-c(2,3)
      c(3,3)=-c(3,3)
      return
   20 c(1,1)=-c(1,1)
      c(2,1)=-c(2,1)
      c(3,1)=-c(3,1)
      c(1,3)=-c(1,3)
      c(2,3)=-c(2,3)
      c(3,3)=-c(3,3)
      return
   30 c(1,1)=-c(1,1)
      c(2,1)=-c(2,1)
      c(3,1)=-c(3,1)
      c(1,2)=-c(1,2)
      c(2,2)=-c(2,2)
      c(3,2)=-c(3,2)
      return
      end
      SUBROUTINE PARSECIF (CIFNAM,NUNIT,CRR)
      DIMENSION CRR(3,3),CELL(7)
      CHARACTER*32 CIFNAM,CR(3,3),cel(7)
      CHARACTER*27 ORIENSTRG11,ORIENSTRG12,ORIENSTRG13,
     1ORIENSTRG21,ORIENSTRG22,ORIENSTRG23,
     2ORIENSTRG31,ORIENSTRG32,ORIENSTRG33,
     3CELLA,CELLB,CELLC,CELLAL,CELLBE,CELLGA,CELLV
        COMMON/DATC/CELL
	CHARACTER*27 TESTSTRG,RELEM
	PARAMETER (ORIENSTRG11='_diffrn_orient_matrix_UB_11')
	PARAMETER (ORIENSTRG12='_diffrn_orient_matrix_UB_12')
	PARAMETER (ORIENSTRG13='_diffrn_orient_matrix_UB_13')
	PARAMETER (ORIENSTRG21='_diffrn_orient_matrix_UB_21')
	PARAMETER (ORIENSTRG22='_diffrn_orient_matrix_UB_22')
	PARAMETER (ORIENSTRG23='_diffrn_orient_matrix_UB_23')
	PARAMETER (ORIENSTRG31='_diffrn_orient_matrix_UB_31')
	PARAMETER (ORIENSTRG32='_diffrn_orient_matrix_UB_32')
	PARAMETER (ORIENSTRG33='_diffrn_orient_matrix_UB_33')
	PARAMETER (CELLA='_cell_length_a')
	PARAMETER (CELLB='_cell_length_b')
	PARAMETER (CELLC='_cell_length_c')
	PARAMETER (CELLAL='_cell_angle_alpha')
	PARAMETER (CELLBE='_cell_angle_beta ')
	PARAMETER (CELLGA='_cell_angle_gamma')
        PARAMETER (CELLV='_cell_volume')
C 
  5	READ(NUNIT,105, END=6) TESTSTRG,RELEM
	IF (TESTSTRG.EQ.ORIENSTRG11) CR(1,1)=RELEM
	IF (TESTSTRG.EQ.ORIENSTRG12) CR(1,2)=RELEM
	IF (TESTSTRG.EQ.ORIENSTRG13) CR(1,3)=RELEM	
	IF (TESTSTRG.EQ.ORIENSTRG21) CR(2,1)=RELEM
	IF (TESTSTRG.EQ.ORIENSTRG22) CR(2,2)=RELEM
	IF (TESTSTRG.EQ.ORIENSTRG23) CR(2,3)=RELEM
	IF (TESTSTRG.EQ.ORIENSTRG31) CR(3,1)=RELEM		
	IF (TESTSTRG.EQ.ORIENSTRG32) CR(3,2)=RELEM
	IF (TESTSTRG.EQ.ORIENSTRG33) CR(3,3)=RELEM
	IF (TESTSTRG.EQ.CELLA)  CEL(1)=RELEM
	IF (TESTSTRG.EQ.CELLB)  CEL(2)=RELEM
	IF (TESTSTRG.EQ.CELLC)  CEL(3)=RELEM
        GO TO 5
  6     rewind nunit
  7     READ(NUNIT,106, END=8) TESTSTRG,RELEM
        IF (TESTSTRG.EQ.CELLA)  CEL(1)=RELEM
	IF (TESTSTRG.EQ.CELLB)  CEL(2)=RELEM
	IF (TESTSTRG.EQ.CELLC)  CEL(3)=RELEM
        go to 7
  8     rewind nunit
  9     READ(NUNIT,107, END=10) TESTSTRG,RELEM
        IF (TESTSTRG.EQ.CELLAL) CEL(4)=RELEM
	IF (TESTSTRG.EQ.CELLBE) CEL(5)=RELEM
	IF (TESTSTRG.EQ.CELLGA) CEL(6)=RELEM
        go to 9
 10     rewind nunit
 11     READ(NUNIT,108, END=12) TESTSTRG,RELEM
        IF (TESTSTRG.EQ.CELLV)  CEL(7)=RELEM
        go to 11
 12     do 13 i=1,3
        do 13 j=1,3
 13     read(cr(i,j),*) crr(i,j)
        do 14 k=1,7
 14     read(cel(k),*) cell(k)
C
        
C
105   FORMAT (A27, A27)
106   FORMAT (A14, A27)
107   FORMAT (A17, A27)
108   FORMAT (A12, A27)
110   FORMAT (3(3A27/))
111   FORMAT (9A27)
115   FORMAT (9F27.10)
120   FORMAT (3(3F20.9/))
200   RETURN
      END
      SUBROUTINE PARSEP4P(NUNIT,CELL)
      DIMENSION ACELL(7),CELL(7)
      CHARACTER*12 ACELL
      CHARACTER*5 S1,S2,S3,S4,S5,S6,S7,S8,CELTEST(2)
 10   READ(NUNIT,100) CELTEST,ACELL
100   FORMAT(A4,A2,6A10,A11,/)
      IF(CELTEST(1).NE.'CELL') GO TO 10
        do 14 k=1,7
 14     read(acell(k),*) cell(k)
      RETURN
      END
