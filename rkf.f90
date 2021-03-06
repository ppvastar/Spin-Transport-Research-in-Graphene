!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     C
!C     ADAPTIVE STEPSIZE DRIVER FOR ODE AND INTEGRATION            C
!C     BY RUNGE-KUTTA-FEHLBERG METHOD                      C
!C__________________________________________________________________________C
!C     ref.Z pp522-532 93/1/12/ 5:15 pm.   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE RKF45(F,NEQ,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK,IWORK)
      integer NEQ
      external F
      real*8 Y(NEQ),T, TOUT, RELERR, ABSERR, WORK(6*NEQ+3)
      integer IFLAG, IWORK(5)
      integer K1M, K1, K2, K3, K4, K5, K6
      K1M=NEQ+1
      K1=K1M+1
      K2=K1+NEQ
      K3=K2+NEQ
      K4=K3+NEQ
      K5=K4+NEQ
      K6=K5+NEQ
      CALL RKFS(F,NEQ,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK(1),&
          WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),&
          WORK(K6),WORK(K6+1),IWORK(1),IWORK(2),IWORK(3),IWORK(4),&
          IWORK(5))
      RETURN
END SUBROUTINE RKF45

SUBROUTINE RKFS(F,NEQ,Y,T,TOUT,RELERR,ABSERR,IFLAG,YP,H,&
    F1,F2,F3,F4,F5,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG,KFLAG)
      external F
      integer NEQ
      real*8 Y(NEQ),YP(NEQ),F1(NEQ),F2(NEQ),F3(NEQ),F4(NEQ), F5(NEQ)
      real*8 T, TOUT, RELERR, ABSERR, H
      real*8 SAVRE, SAVAE
      integer NFE, KOP, INIT, JFLAG, KFLAG, IFLAG
      LOGICAL HFAILD,OUTPUT
      integer ISTEP
      real*8 EPS, EPSP1, U26, RER, A, DT, TOLN, TOL, YPK
      real*8 SCALE, AE, HMIN, EEOET, ET, EE, ESTTOL, S
      COMMON/ISTEP/ISTEP
      real*8 REMIN
      DATA REMIN/1.D-12/
      integer MAXNFE
      DATA MAXNFE/3000/
      integer MFLAG, K
      
      IF(NEQ.LT.1) GOTO 10
      IF((RELERR.LT.0.0).OR.(ABSERR.LT.0.0)) GOTO 10
      MFLAG=IABS(IFLAG)
      IF((MFLAG.EQ.0).OR.(MFLAG.GT.8)) GOTO 10
      IF(MFLAG.NE.1) GOTO 20
      EPS=1.D0
 5    EPS=EPS*0.5D0
      EPSP1=EPS+1.D0
      IF(EPSP1.GT.1.D0) GOTO 5
      U26=26.D0*EPS
      GOTO 50
 10   IFLAG=8
      RETURN
 20   IF((T.EQ.TOUT).AND.(KFLAG.NE.3)) GOTO 10
      IF(MFLAG.NE.2) GOTO 25
      IF((KFLAG.EQ.3).OR.(INIT.EQ.0)) GOTO 45
      IF(KFLAG.EQ.4) GOTO 40
      IF((KFLAG.EQ.5).AND.(ABSERR.EQ.0.D0)) GOTO 30
      IF((KFLAG.EQ.6).AND.(RELERR.LE.SAVRE).AND.(ABSERR.LE.SAVAE)) GOTO 30
      GOTO 50
 25   IF(IFLAG.EQ.3) GOTO 45
      IF(IFLAG.EQ.4) GOTO 40
      IF((IFLAG.EQ.5).AND.(ABSERR.GT.0.D0)) GOTO 45
 30   STOP 12
 40   NFE=0
      IF(MFLAG.EQ.2) GOTO 50
 45   IFLAG=JFLAG
      IF(KFLAG.EQ.3) MFLAG=IABS(IFLAG)
 50   JFLAG=IFLAG
      KFLAG=0
      SAVRE=RELERR
      SAVAE=ABSERR
      RER=2.D0*EPS+REMIN
      IF(RELERR.GE.RER) GOTO 55
      RELERR=RER
      IFLAG=3
      KFLAG=3
      RETURN
 55   DT=TOUT-T
      IF(MFLAG.EQ.1) GOTO 60
      IF(INIT.EQ.0) GOTO 65
      GOTO 80
 60   INIT=0
      KOP=0
      A=T
      CALL F(A,Y,YP)
      NFE=1
      IF(T.NE.TOUT) GOTO 65
      IFLAG=2
      RETURN
 65   INIT=1
      H=DABS(DT)
      TOLN=0.D0
      DO 70 K=1,NEQ
         TOL=RELERR*DABS(Y(K))+ABSERR
         IF(TOL.LE.0.D0) GOTO 70
         TOLN=TOL
         YPK=DABS(YP(K))
         IF(YPK*H**5.GT.TOL) H=(TOL/YPK)**0.2D0
 70   CONTINUE
      IF(TOLN.LE.0.D0) H=0.D0
      H=DMAX1(H,U26*DMAX1(DABS(T),DABS(DT)))
      JFLAG=ISIGN(2,IFLAG)
 80   H=DSIGN(H,DT)
      IF(DABS(H).GE.2.D0*DABS(DT)) KOP=KOP+1
      IF(KOP.NE.100) GOTO 85
      KOP=0
      IFLAG=7
      RETURN
 85   IF(DABS(DT).GT.U26*DABS(T)) GOTO 95
      DO 90 K=1,NEQ
         Y(K)=Y(K)+DT*YP(K)
 90   END DO
      A=TOUT
      CALL F(A,Y,YP)
      NFE=NFE+1
      GOTO 300
 95   OUTPUT=.FALSE.
      SCALE=2.D0/RELERR
      AE=SCALE*ABSERR
 100  HFAILD=.FALSE.
      HMIN=U26*DABS(T)
      DT=TOUT-T
      IF(DABS(DT).GE.2.D0*DABS(H)) GOTO 200
      IF(DABS(DT).GT.DABS(H)) GOTO 150
      OUTPUT=.TRUE.
      H=DT
      GOTO 200
 150  H=0.5D0*DT
 200  IF(NFE.LE.MAXNFE) GOTO 220
      IFLAG=4
      KFLAG=4
      RETURN
 220  CALL FEHL(F,NEQ,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE=NFE+5
      EEOET=0.D0
      DO 250 K=1,NEQ
         ET=DABS(Y(K))+DABS(F1(K))+AE
         IF(ET.GT.0.D0) GOTO 240
         IFLAG=5
         RETURN
240      EE=DABS((-2090.D0*YP(K)+(21970.D0*F3(K)-15048.D0*F4(K)))&
             +(22528.D0*F2(K)-27360.D0*F5(K)))
         EEOET=DMAX1(EEOET,EE/ET)
 250  END DO
      ESTTOL=DABS(H)*EEOET*SCALE/752400.D0
      IF(ESTTOL.LE.1.D0) GOTO 260
      HFAILD=.TRUE.
      OUTPUT=.FALSE.
      S=0.1D0
      IF(ESTTOL.LT.59049.D0) S=0.9D0/ESTTOL**0.2D0
      H=S*H
      IF(DABS(H).GT.HMIN) GOTO 200
      IFLAG=6
      KFLAG=6
      RETURN
 260  T=T+H
      ISTEP=ISTEP+1
      DO 270 K=1,NEQ
         Y(K)=F1(K)
 270  END DO
      A=T
      CALL F(A,Y,YP)
      NFE=NFE+1
      S=5.D0
      IF(ESTTOL.GT.1.889568D-4) S=0.9D0/ESTTOL**0.2D0
      IF(HFAILD) S=DMIN1(S,1.D0)
      H=DSIGN(DMAX1(S*DABS(H),HMIN),H)
      IF(OUTPUT) GOTO 300
      IF(IFLAG.GT.0) GOTO 100
      IFLAG=-2
      RETURN
 300  T=TOUT
      IFLAG=2
      RETURN
     END 

     SUBROUTINE FEHL(F,NEQ,Y,T,H,YP,F1,F2,F3,F4,F5,S)
!c$$$  IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F
      integer NEQ
      real*8 Y(NEQ),YP(NEQ),F1(NEQ),F2(NEQ),F3(NEQ),F4(NEQ),F5(NEQ),S(NEQ), H, T
      real*8 CH
      integer K
      CH=H/4.D0
      DO 221 K=1,NEQ
         F5(K)=Y(K)+CH*YP(K)
 221  END DO
      CALL F(T+CH,F5,F1)
      CH=3.D0*H/32.D0
      DO 222 K=1,NEQ
         F5(K)=Y(K)+CH*(YP(K)+3.D0*F1(K))
 222  END DO
      CALL F(T+3.D0*H/8.D0,F5,F2)
      CH=H/2197.D0
      DO 223 K=1,NEQ
         F5(K)=Y(K)+CH*(1932.D0*YP(K)+(7296.D0*F2(K)-7200.D0*F1(K)))
 223  END DO
      CALL F(T+12.D0*H/13.D0,F5,F3)
      CH=H/4104.D0
      DO 224 K=1,NEQ
         F5(K)=Y(K)+CH*((8341.D0*YP(K)-845.D0*F3(K))+(29440.D0*F2(K)-32832.D0*F1(K)))
 224 END DO
      CALL F(T+H,F5,F4)
      CH=H/20520.D0
      DO 225 K=1,NEQ
         F1(K)=Y(K)+CH*((-6068.D0*YP(K)+(9295.D0*F3(K)-5643.D0*F4(K)))&
             +(41040.D0*F1(K)-28352.D0*F2(K)))
 225  END DO
      CALL F(T+H/2.D0,F1,F5)
      CH=H/7618050.D0
      DO 230 K=1,NEQ
       S(K)=Y(K)+CH*((902880.D0*YP(K)+(3855735.D0*F3(K)-1371249.D0&
           *F4(K)))+(3953664.D0*F2(K)+277020.D0*F5(K)))
 230  END DO
      RETURN
END SUBROUTINE FEHL


