!  LAPJVTEST.f90 
!
!  FUNCTIONS:
!  LAPJVTEST - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: LAPJVTEST
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program LAPJVTEST
    use global
    implicit none

    ! Variables
    character(150) :: fileName !job file name
    integer :: n, i = 1, j = 1   
    integer, dimension(100) :: row
    integer, dimension(100,100) :: c
    integer :: x(100), &                !col assigned to row
               y(100), &                !row assigned to column
               u(100), &                !Dual Row variable
               v(100), &                !Dual column variable
               z 
    ! Body of LAPJVTEST
    
    !!get the file name from program argument
    call getarg(1, fileName)
    fileName = trim(fileName)
    
    !read the jobfile
    call readjobfile(fileName)
    
    !print the jobtable
    do i = 1, numjob
        write(*,*) jobtable(i).jobi, jobtable(i).pi, jobtable(i).di, jobtable(i).wi, jobtable(i).ri   
    end do
    
    !read the cost matrix file into cost matrix and get total number of jobs
    !call readcostmatrix(fileName, n, c)
    
    !print the cost matrix
    !j = 1
    !do i=1, n
    !    write(*, '(4I5)') (c(i,j), j=1,n)   
    !end do
    
    !run LAPJV
    !call JOVOFDTEST(n, c, x, y, u, v, z)
    
    !print result
    !write(*, '(A15, I5)')  "Optimal value is ", z
        
    write(*,*) "Please press Enter to continue with next case......"
    read(*,*)   

end program LAPJVTEST
    
