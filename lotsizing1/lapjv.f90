module lapjv

    INTEGER :: maxDim = 100
    
contains    
subroutine jovofd(N,C,X,Y,U,V,Z)
    implicit none
      INTEGER C(100,100),X(100),Y(100),U(100),V(100)
      INTEGER H,Z,L0,V0,VJ,DJ,UP,LOW
      INTEGER LAB(100),D(100),FREE(100),COL(100)
    INTEGER N,MIN,LAST,K,J1,J0,I,CNT,L,J,I0

! INITIALIZATION
      DO 10 I=1,N
        X(I)=0
   10 CONTINUE
      DO 20 J0=1,N
        J=N-J0+1
        VJ=C(J,1)
        I0=1
        DO 15 I=2,N
          IF (C(J,I).LT.VJ) THEN
            VJ=C(J,I)
            I0=I
          END IF
   15   CONTINUE
        V(J)=VJ
        COL(J)=J
        IF (X(I0).EQ.0) THEN
          X(I0)=J
          Y(J)=I0
        ELSE
          X(I0)=-ABS(X(I0))
          Y(J)=0
        END IF
   20 CONTINUE
      L=0
      DO 40 I=1,N
        IF (X(I).EQ.0) THEN
          L=L+1
          FREE(L)=I
          GOTO 40
        END IF
        IF (X(I).LT.0) THEN
          X(I)=-X(I)
        ELSE
          J1=X(I)
          MIN=1.E14
          DO 31 J=1,N
            IF (J.EQ.J1) GOTO 31
            IF (C(J,I)-V(J).LT.MIN) MIN=C(J,I)-V(J)
   31     CONTINUE
          V(J1)=V(J1)-MIN
        END IF
   40 CONTINUE
! IMPROVE THE INITIAL SOLUTION
      CNT=0
      IF (L.EQ.0) GOTO 1000
   41 L0=L
      K=1
      L=0
   50 I=FREE(K)
      K=K+1
      V0=C(1,I)-V(1)
      J0=1
      VJ=1.E14
      DO 42 J=2,N
        H=C(J,I)-V(J)
        IF (H.LT.VJ) THEN
        IF (H.GE.V0) THEN
          VJ=H
          J1=J
        ELSE
          VJ=V0
          V0=H
          J1=J0
          J0=J
        END IF
      END IF
   42 CONTINUE
      I0=Y(J0)
      IF (V0.LT.VJ) THEN
        V(J0)=V(J0)-VJ+V0
      ELSE
        IF (I0.EQ.0) GOTO 43
        J0=J1
        I0=Y(J1)
      END IF
      IF (I0.EQ.0) GOTO 43
      IF (V0.LT.VJ) THEN
        K=K-1
        FREE(K)=I0
      ELSE
        L=L+1
        FREE(L)=I0
      END IF
   43 X(I)=J0
      Y(J0)=I
      IF (K.LE.L0) GOTO 50
      CNT=CNT+1
      IF ((L.GT.0).AND.(CNT.LT.2)) GOTO 41
! AUGMENTATION PART
      L0=L
      DO 90 L=1,L0
        I0=FREE(L)
        DO 51 J=1,N
          D(J)=C(J,I0)-V(J)
          LAB(J)=I0
   51   CONTINUE
        UP=1
        LOW=1
   60   LAST=LOW-1
        MIN=D(COL(UP))
        UP=UP+1
        DO 52 K=UP,N
          J=COL(K)
          DJ=D(J)
          IF (DJ.LE.MIN) THEN
            IF (DJ.LT.MIN) THEN
              MIN=DJ
              UP=LOW
            END IF
            COL(K)=COL(UP)
            COL(UP)=J
            UP=UP+1
        END IF
   52   CONTINUE
        DO 53 H=LOW,UP-1
          J=COL(H)
          IF (Y(J).EQ.0) GOTO 70
   53   CONTINUE
   55   J0=COL(LOW)
        LOW=LOW+1
        I=Y(J0)
        H=C(J0,I)-V(J0)-MIN
        DO 62 K=UP,N
          J=COL(K)
          VJ=C(J,I)-V(J)-H
          IF (VJ.LT.D(J)) THEN
            D(J)=VJ
            LAB(J)=I
            IF (VJ.EQ.MIN) THEN
              IF (Y(J).EQ.0) GOTO 70
              COL(K)=COL(UP)
              COL(UP)=J
              UP=UP+1
            END IF
          END IF
   62   CONTINUE
        IF (LOW.EQ.UP) GOTO 60
        GOTO 55
   70   DO 71 K=1,LAST
          J0=COL(K)
          V(J0)=V(J0)+D(J0)-MIN
   71   CONTINUE
   80   I=LAB(J)
        Y(J)=I
        K=J
        J=X(I)
        X(I)=K
        IF (I0.NE.I) GOTO 80
   90 CONTINUE
      Z=0
      DO 100 I=1,N
        U(I)=C(X(I),I)-V(X(I))
        Z=Z+C(X(I),I)
  100 CONTINUE
 1000 RETURN

end subroutine jovofd 

SUBROUTINE JOVOSAP(N,CC,KK,FIRST,X,Y,U,V,H)
      INTEGER CC(2500),KK(2500),FIRST(101),X(100),Y(100),U(100),V(100)
      INTEGER H,CNT,L0,T,T0,TD,V0,VJ,DJ
      INTEGER LAB(100),D(100),FREE(100),TODO(100)
      LOGICAL OK(100)

!C THIS SUBROUTINE SOLVES THE SPARSE LINEAR ASSIGNMENT PROBLEM
!C ACCORDING TO R. JONKER & A. VOLGENANT,
!C UNIVERSITY OF AMSTERDAM,
!C COMPUTING 38, 325-340 (1987),
!C A SHORTEST AUGMENTING PATH ALGORITHM
!C FOR DENSE AND SPARSE LINEAR ASSIGNMENT PROBLEMS.
!C
!C INPUT PARAMETERS :
!C N = NUMBER OF ROWS AND COLUMNS
!C C = WEIGHT MATRIX
!C
!C OUTPUT PARAMETERS
!C X = COL ASSIGNED TO ROW
!C Y = ROW ASSIGNED TO COL
!C U = DUAL ROW VARIABLE
!C V = DUAL COLUMN VARIABLE
!C H = VALUE OF OPTIMAL SOLUTION
!C
!C INITIALIZATION
      DO 10 J=1,N
      V(J)=1.E14
   10 CONTINUE
      DO 20 I=1,N
      X(I)=0
      DO 15 T=FIRST(I),FIRST(I+1)-1
      J=KK(T)
      IF (CC(T).LT.V(J)) THEN
        V(J)=CC(T)
        Y(J)=I
      END IF
   15 CONTINUE
   20 CONTINUE
      DO 30 J=1,N
      J0=N-J+1
      I=Y(J0)
      IF (X(I).NE.0) THEN
        X(I)=-ABS(X(I))
        Y(J0)=0
      ELSE
        X(I)=J0
      END IF
   30 CONTINUE
      L=0
      DO 40 I=1,N
      IF (X(I).EQ.0) THEN
        L=L+1
        FREE(L)=I
        GOTO 40
      END IF
      IF (X(I).LT.0) THEN
        X(I)=-X(I)
      ELSE
        J1=X(I)
        MIN=1.E14
        DO 31 T=FIRST(I),FIRST(I+1)-1
          J=KK(T)
          IF (J.EQ.J1) GOTO 31
          IF (CC(T)-V(J).LT.MIN) MIN=CC(T)-V(J)
   31   CONTINUE
        V(J1)=V(J1)-MIN
      END IF
   40 CONTINUE
! IMPROVE THE INITIAL SOLUTION
      CNT=0
      IF (L.EQ.0) GOTO 1000
   41 L0=L
      K=1
      L=0
   50 I=FREE(K)
      K=K+1
      V0=1.E14
      VJ=1.E14
      DO 42 T=FIRST(I),FIRST(I+1)-1
      J=KK(T)
      H=CC(T)-V(J)
      IF (H.LT.VJ) THEN
        IF (H.GE.V0) THEN
          VJ=H
          J1=J
        ELSE
          VJ=V0
          V0=H
          J1=J0
          J0=J
        END IF
      END IF
   42 CONTINUE
      I0=Y(J0)
      IF (V0.LT.VJ) THEN
        V(J0)=V(J0)-VJ+V0
      ELSE
        IF (I0.EQ.0) GOTO 43
        J0=J1
        I0=Y(J1)
      END IF
      IF (I0.EQ.0) GOTO 43
      IF (V0.LT.VJ) THEN
        K=K-1
        FREE(K)=I0
      ELSE
        L=L+1
        FREE(L)=I0
      END IF
   43 X(I)=J0
      Y(J0)=I
      IF (K.LE.L0) GOTO 50
      CNT=CNT+1
      IF ((L.GT.0).AND.(CNT.LT.2)) GOTO 41
! AUGMENTATION PART
      L0=L
      DO 90 L=1,L0
      DO 51 J=1,N
      OK(J)=.FALSE.
      D(J)=1.E14
   51 CONTINUE
      MIN=1.E14
      I0=FREE(L)
      TD=N
      DO 52 T=FIRST(I0),FIRST(I0+1)-1
      J=KK(T)
      DJ=CC(T)-V(J)
      D(J)=DJ
      LAB(J)=I0
      IF (DJ.LE.MIN) THEN
        IF (DJ.LT.MIN) THEN
          MIN=DJ
          K=1
          TODO(1)=J
        ELSE
          K=K+1
          TODO(K)=J
        END IF
      END IF
   52 CONTINUE
      DO 53 H=1,K
      J=TODO(H)
      IF (Y(J).EQ.0) GOTO 80
      OK(J)=.TRUE.
   53 CONTINUE
! REPEAT UNTIL A FREE ROW HAS BEEN FOUND
   60 J0=TODO(K)
      K=K-1
      I=Y(J0)
      TODO(TD)=J0
      TD=TD-1
      T0=FIRST(I)
      T=T0-1
   61 T=T+1
      IF (KK(T).NE.J0) GOTO 61
      H=CC(T)-V(J0)-MIN
      DO 62 T=T0,FIRST(I+1)-1
      J=KK(T)
      IF (.NOT. OK(J)) THEN
        VJ=CC(T)-H-V(J)
        IF (VJ.LT.D(J)) THEN
          D(J)=VJ
          LAB(J)=I
          IF (VJ.EQ.MIN) THEN
            IF (Y(J).EQ.0) GOTO 70
            K=K+1
            TODO(K)=J
            OK(J)=.TRUE.
          END IF
        END IF
      END IF
   62 CONTINUE
      IF (K.NE.0) GOTO 60
      MIN=1.E14-1
      DO 63 J=1,N
      IF (D(J).LE.MIN) THEN
        IF (.NOT. OK(J)) THEN
          IF (D(J).LT.MIN) THEN
            MIN=D(J)
            K=1
            TODO(1)=J
          ELSE
            K=K+1
            TODO(K)=J
          END IF
        END IF
      END IF
   63 CONTINUE
      DO 64 J0=1,K
      J=TODO(J0)
      IF (Y(J).EQ.0) GOTO 70
      OK(J)=.TRUE.
   64 CONTINUE
      GOTO 60
   70 IF (MIN.EQ.0) GOTO 80
      DO 71 K=TD+1,N
      J0=TODO(K)
      V(J0)=V(J0)+D(J0)-MIN
   71 CONTINUE
   80 I=LAB(J)
      Y(J)=I
      K=J
      J=X(I)
      X(I)=K
      IF (I0.NE.I) GOTO 80
   90 CONTINUE
      H=0
      DO 100 I=1,N
      J=X(I)
      T=FIRST(I)-1
  101 T=T+1
      IF (KK(T).NE.J) GOTO 101
      DJ=CC(T)
      U(I)=DJ-V(J)
      H=H+DJ
  100 CONTINUE
1000 return
end subroutine JOVOSAP


    end module 