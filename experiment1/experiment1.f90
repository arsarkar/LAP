!  experiment1.f90 
!
!  FUNCTIONS:
!  experiment1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: experiment1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program experiment1
    use casegenerator
    use global
    use BranchBound
    implicit none
    double precision :: r=0.6, t=0.6
    !!open the file to write
    open(unit = output, file = "output.txt")
    
    !!configure the data generator
    call configure(5,10,1,3)
    numjob = 5 
    
    !!generate the data
    call potts1982(numjob,r,t)
    
    !!get jobs
    allocate(jobtable(numjob))
    jobtable(1:numjob) = jobs(1:numjob)
    
    !print jobs
    write(output,'(A17,A150)') "jobs created for r=0.6, t=0.6"
    call printjobs()
    
    !read the cost matrix file into cost matrix and get total number of jobs
    call assignmatrix()
    
    !print cost matrix
    write(output,'(A21)') "cost matrix assigned " 
    !call printcostmatrix()
    
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
    
    end program experiment1

