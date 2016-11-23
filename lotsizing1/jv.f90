!  jv.f90 
!
!  FUNCTIONS:
!  jv - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: jv
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    module jv
        implicit none
        !character(150) :: fileName, outputfile !job file name
        !integer :: n, i = 1, j = 1   
        !integer, dimension(100) :: row
        !integer :: x(100), &                !col assigned to row
        !           y(100), &                !row assigned to column
        !           u(100), &                !Dual Row variable
        !           v(100), &                !Dual column variable
        !           z 
        !Integer :: c(100,100)
        !integer :: output = 5
        integer :: maxdim = 100
    
    contains
    
    subroutine solveLAP(costmat,r2c,c2r,cost)
        use global
        implicit none
        
        integer :: costmat(dim,dim),c2r(dim),r2c(dim)
        integer, dimension(:,:), allocatable :: CM
        integer, dimension(:), allocatable :: X,Y,U,V
        integer :: cost, allocateStat
        
        !assign the maximum size for JOVOFD 
        do while (dim > maxdim)
            maxdim = maxdim + 100
        end do
        
        !reallocate CM
        allocate(CM(maxdim, maxdim), STAT = allocateStat)
        allocate(X(maxdim), STAT = allocateStat)
        allocate(Y(maxdim), STAT = allocateStat)
        allocate(U(maxdim), STAT = allocateStat)
        allocate(V(maxdim), STAT = allocateStat)
        
        !copy costmatrix
        CM(1:dim,1:dim) = costmat(1:dim,1:dim)
        
        !solve the LAP
        call jovofd(dim,CM,X,Y,U,V,cost)
        
        !copy back the solutions
        r2c(1:dim) = Y(1:dim)
        c2r(1:dim) = X(1:dim)        
    
    end subroutine solveLAP
        
    !****************************************************************************
    !
    !  SUBROUTINE: JOVOFD
    !
    !  PURPOSE:  This Subroutine Solves The Full Density Linear Assignment Problem
    !            According To "A Shortest Augmenting Path Algorithm For Dense And 
    !            Sparse Linear Assignment Problems," Computing 38, 325-340, 1987
    !            By R. Jonker And A. Volgenant, University Of Amsterdam.
    !  
    !****************************************************************************
    !  INPUT: N = Number Of Rows And Columns
    !         C = Weight Matrix
    !
    !  Output: X = Col Assigned To Row
    !          Y = Row Assigned To Col
    !          U = Dual Row Variable
    !          V = Dual Column Variable
    !          Z = Value Of Optimal Solution
    !****************************************************************************   
    !=================================================
    !print solved matrix 
    !=================================================
    subroutine jovofd(N,C,X,Y,U,V,Z)
        implicit none
          INTEGER C(maxdim,maxdim),X(maxdim),Y(maxdim),U(maxdim),V(maxdim)
          INTEGER H,Z,L0,V0,VJ,DJ,UP,LOW
          INTEGER LAB(maxdim),D(maxdim),FREE(maxdim),COL(maxdim)
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
    
    end module
