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
    use BranchBound
    use jv
    implicit none

    ! Variables
    character(150) :: fileName, outputfile !job file name
    integer :: n, i = 1, j = 1   
    integer, dimension(100) :: row
    integer :: x(100), &                !col assigned to row
               y(100), &                !row assigned to column
               u(100), &                !Dual Row variable
               v(100), &                !Dual column variable
               z 
    ! Body of LAPJVTEST
    !!get the file name from program argument
    call getarg(1, fileName)
    fileName = trim(fileName)
    
    !!get output filename 
    call getarg(2, outputfile)
    outputfile = trim(outputfile)
    
    !!open the file to write
    open(unit = output, file = outputfile)
    
    !read the jobfile
    write(output,'(A13,A150)') "Reading file ", filename 
    write(*,'(A13,A150)') "Reading file ", filename
    call readjobfile(fileName)
    
    !print jobs
    write(output,'(A17,A150)') "jobs parsed from ", filename 
    write(*,'(A17,A150)') "jobs parsed from ", filename 
    call printjobs()
    
    !read the cost matrix file into cost matrix and get total number of jobs
    call assignmatrix()
    
    !print cost matrix
    write(output,'(A21)') "cost matrix assigned " 
    write(*,'(A21)') "cost matrix assigned " 
    !call printcostmatrix()
    
    allocate(heuristicSchedule(dim))
    
    !apply heuristic to find upper bound
    call wsrpt(upperBound, heuristicSchedule)
    
    !test branch and bound
    call DFSBB(c)
    !call printcostmatrix(c)
    
    !run LAPJV
    !call JOVOFDTEST(n, c%, x, y, u, v, z)
    
    !print result
    !write(*, '(A15, I5)')  "Optimal value is ", z
    !call printsolvedmatrix(c, x, y)
        
    write(*,*) "Please press Enter to continue with next case......"
    read(*,*)   
    
end program LAPJVTEST
    
