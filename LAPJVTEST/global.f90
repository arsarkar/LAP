!A global module containing types and functions 
module global

    type jobstruct
        integer :: jobi
        integer :: pi
        integer :: di
        integer :: wi
        integer :: ri
    end type jobstruct    

    type (jobstruct), dimension(:), allocatable :: jobtable
    
    integer :: numjob = 0

contains
!****************************************************************************
!
!  SUBROUTINE: readcostmatrix
!
!  PURPOSE:  This subroutine reads the cost matrix file and returns the number
!            of jobs and the cost matrix
!  
!****************************************************************************
!  INPUT: costfile = cost file path full
!        
!
!  Output:  costmatrix = cost matrix
!           n = number of jobs
!****************************************************************************    
subroutine readcostmatrix(costfile, n, costmatrix)

   character(150) :: costfile, inputline
   integer, dimension(100,100) :: costmatrix
   integer, dimension(100) :: row
   integer :: fileUnit = 15, ioresult, n, i = 1, j = 1
   
    !open the file    
    open(fileUnit,file=costfile, status='old', position='rewind') 
    
    !read the first line
    read(fileUnit, '(A150)', iostat=ioResult) inputline
    inputline = trim(inputline)
    read(inputline, *) n

    
    !read the matrix
    do while (i <= n)
        
        read(fileUnit, '(A150)', iostat=ioResult) inputline
        read(inputline, *) row(1:n)
        j = 1
        do while (j<=n)
           costmatrix(i,j) = row(j)
           j = j + 1
        end do
        
        i = i + 1
        
    end do    

end subroutine readcostmatrix 

!****************************************************************************
!
!  SUBROUTINE: readjobfile
!
!  PURPOSE:  This subroutine reads the job file into jobtable 
!            global jobtable is populated with a size of number of job
!  
!****************************************************************************
!  INPUT: jobfile = full path of job file
!**************************************************************************** 
subroutine readjobfile(jobFile)

    character(150) :: jobfile, inputline
    integer :: fileUnit = 10, ioresult, i = 1, j = 1
    integer, dimension(10) :: row
    type (jobstruct) :: job
    
    !open the file    
    open(fileUnit,file=jobfile, status='old', position='rewind') 
    
    !read the first line
    read(fileUnit, '(A150)', iostat=ioResult) inputline
    !trim(inputline)
    read(inputline, *) numjob
    
    !redim the array of jobs to number of jobs
    allocate(jobtable(numjob))
    
    !read the job matrix
    do while (i <= numjob)
        read(fileUnit, '(A150)', iostat=ioResult) inputline
        read(inputline, *) job
        jobtable(i) = job
        i = i + 1
    end do    

end subroutine readjobfile
    
subroutine matrixtransfer(jobfile, n, c)

    character(150) :: jobfile
    integer :: fileUnit = 10, ioresult, n, i = 1, j = 1, p = 0, d = 0, dim = 0
    integer, dimension(100) :: job
    integer, dimension(100,100) :: c
    !open the file    
    open(fileUnit,file=jobfile, status='old', position='rewind') 
    
    !read the first line
    read(fileUnit, '(A150)', iostat=ioResult) inputline
    !trim(inputline)
    read(inputline, *) n
    
    !calculate the total dimensation fo cost matrix
    ! summation of all processing time as time is granulated as unit time
    do while (j <= n)
        read(fileUnit, '(A150)', iostat=ioResult) inputline
        read(inputline, *) job(1:5)
        dim = dim + job(2)
    end do    
    
    do while (j <= n)  
        
        read(fileUnit, '(A150)', iostat=ioResult) inputline
        read(inputline, *) job(1:5)
        p = 0
        d = 0
        
        !while (p <= row(2)) !add equal number of rows for each process 
            
            
        
        !end while
        
    end do  
    

end subroutine matrixtransfer
   
    
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
subroutine JOVOFDTEST(N,C,X,Y,U,V,Z)
    
    implicit none
    
    !variable declaration
    integer :: C(100,100), &            !cost/weight matrix
               X(100), &                !col assigned to row
               Y(100), &                !row assigned to column
               U(100), &                !Dual Row variable
               V(100), &                !Dual column variable
               N, &                     !number of row and column
               H, &                     !
               Z, &                     !value of optimal solution
               L0, &                    
               V0, &                    
               VJ, &                    
               DJ, &                    
               UP, &                   
               LOW, &
               MIN, &
               LAB(100), &
               D(100), &
               FREE(100), &
               COL(100)
    
    !counter variables
    integer:: I, I0, J, J0, J1, K, L, CNT, LAST

    !initialization
    
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
      IF (L.EQ.0) return
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

end subroutine JOVOFDTEST    
    
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
subroutine JOVOFD(N,C,X,Y,U,V,Z)
    
    implicit none
    
    !variable declaration
    integer :: C(100,100), &            !cost/weight matrix
               X(100), &                !col assigned to row
               Y(100), &                !row assigned to column
               U(100), &                !Dual Row variable
               V(100), &                !Dual column variable
               N, &                     !number of row and column
               H, &                     !
               Z, &                     ! value of optimal solution
               L0, &                    
               V0, &                    
               VJ, &                    
               DJ, &                    
               UP, &                   
               LOW, &
               MIN, &
               LAB(100), &
               D(100), &
               FREE(100), &
               COL(100)
    
    !loop variables
    integer:: I, I0, J, J0, J1, L

    !initialization
    
    !initialize row assignment matrix
    do I=1,N 
       X(I)=0
    end do
    
    do J0=1,N
       J=N-J0+1
       VJ=C(J,1)
       I0=1
       do I=2,N
          if (C(J,I) .lt. VJ) then
             VJ=C(J,I)
             I0=I
          end if
       end do    
       V(J)=VJ
       COL(J)=J
       if (X(I0) .eq. 0) then
          X(I0)=J
          Y(J)=I0
       else
          X(I0)=-ABS(X(I0))
          Y(J)=0
       end if
    end do   
    L=0
    do I=1,N
       if (X(I) .eq. 0) then
           L=L+1
           FREE(L)=I
           continue
       end if
       if (X(I) .lt. 0) then
           X(I)=-X(I)
       else
           J1=X(I)
           MIN=1.E14
           do J=1,N
               if (J .eq. J1) continue
               if (C(J,I)-V(J) .lt. MIN) then
                   MIN=C(J,I)-V(J)
               end if
           end do
           V(J1)=V(J1)-MIN
       end if
    end do
    
   ! !improve the initial solution
   ! CNT=0
   ! IF (L.EQ.0) then
   ! do while ((L.GT.0).AND.(CNT.LT.2))
   !     L0=L
   !     K=1
   !     L=0
   !     do while (K.LE.L0)
   !         I=FREE(K)
   !         K=K+1
   !         V0=C(1,I)-V(1)
   !         J0=1
   !         VJ=1.E14
   !         do J=2,N
   !             H=C(J,I)-V(J)
   !             if (H.LT.VJ) then
   !                 if (H.GE.V0) then
   !                     VJ=H
   !                     J1=J
   !                 else
   !                     VJ=V0
   !                     V0=H
   !                     J1=J0
   !                     J0=J
   !                 end if
   !             end if
   !         end do
   !     I0=Y(J0)
   !     if (V0.LT.VJ) then
   !         V(J0)=V(J0)-VJ+V0
   !     else
   !         if (I0.EQ.0) GOTO 43
   !     J0=J1
   !     I0=Y(J1)
   !   END IF
   !   IF (I0.EQ.0) GOTO 43
   !   IF (V0.LT.VJ) THEN
   !     K=K-1
   !     FREE(K)=I0
   !   ELSE
   !     L=L+1
   !     FREE(L)=I0
   !   END IF
   !43 X(I)=J0
   !   Y(J0)=I
   !   end do
   !   CNT=CNT+1
   ! end do
   !
   ! !augmentation part      
   !   
   ! end if  
end subroutine JOVOFD    
    
end module    