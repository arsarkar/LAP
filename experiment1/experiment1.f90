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
    call configure(maxpi=2,maxwi=10,minwi=1,seed=3)
    numjob = 20
    
    !!generate the data
    call potts1982(numjob,r,t)
    
    !!get jobs
    allocate(jobtable(numjob))
    jobtable(1:numjob) = jobs(1:numjob)
    
    !sort jobs
    call sortjobs()
    
    !print jobs
    write(output,'(A17,A150)') "jobs created for r=0.6, t=0.6"
    call printjobs()
    
    !read the cost matrix file into cost matrix and get total number of jobs
    call assignmatrix()
    
    !print cost matrix
    write(output,'(A21)') "cost matrix assigned " 
    !call printcostmatrix()
    
    allocate(heuristicSchedule(dim))
    
    !apply heuristic to find upper bound
    call wsrpt(upperBound, heuristicSchedule)
    
    !test branch and bound
    call DFSBB(c)
    !call printcostmatrix(c)
        
    write(*,*) "Please press Enter to exit......"
    read(*,*)   
    
    end program experiment1

