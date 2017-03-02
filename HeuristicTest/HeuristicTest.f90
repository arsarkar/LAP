!  HeuristicTest.f90 
!
!  FUNCTIONS:
!  HeuristicTest - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: HeuristicTest
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program HeuristicTest
    use casegenerator
    use global
    use single_premp_p_r_wt    
    implicit none
    
    double precision :: r=0.6, t=0.6
    integer:: i, tt
    
    !open the file to write
    open(unit = output, file = "output.txt")
    
    !configure the data generator 
    call configure(maxpi=2,maxwi=10,minwi=1,seed=3)
    numjob = 15

    !!generate the data
    call potts1982(numjob,r,t)
    
    !!get jobs
    allocate(jobtable(numjob))
    jobtable(1:numjob) = jobs(1:numjob)
    
    !print jobs
    write(output,'(A17,A150)') "jobs created for r=0.6, t=0.6"
    call printjobs()
    
    !call the heuristic
    call slackX1(tt)
    
    !print the schedule
    write(output,'(A25)') "The schedule assigned "
    write(output, 30) (schedule(i), i=1,dim) 
    write(output, '(A20, I5)') "Total tardiness = ", tt
    
    ! Body of HeuristicTest
    print *, 'Hello World'
    
    write(*,*) "Please press Enter to exit......"
    read(*,*) 
    
30  format(<dim>(I3,1X,"|",1X))
    end program HeuristicTest

