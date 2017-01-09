!A global module containing types and functions 
module global

    !type to hold information for one job
    type jobstruct
        integer :: jobi
        integer :: pi
        integer :: di
        integer :: wi
        integer :: ri
    end type jobstruct  
    
    type cmat
        integer, allocatable, dimension(:,:) :: m       !2-D cost matrix
        integer, allocatable, dimension(:) :: r2c       !row to column assignment k-th value is the column assigned to row k
        integer, allocatable, dimension(:) :: c2r       !column to row assignment k-th value is the row assigned to column k
        integer z                                       !cost of this solution
    end type cmat    
    
    !type violation to store one violation
    !ct is the completion time violation and nct is the non completion time violations
    type violation
        integer ct(2)
        integer nct(2)
    end type violation
    
    !type solution
    type solution
        type(cmat) c            !Cost matrix
        integer sr              !solved using strategy 1->including pivot 2-> excluding pivot
        integer pivot(2)        !matrix cell acting as pivot
        integer tolerance       !tolerance of this solution
        integer depth           !depth of this sol in solution tree
    end type solution
    
    !job table
    type (jobstruct), dimension(:), allocatable :: jobtable
    
    !total number of jobs supplied
    integer :: numjob = 0
    
    !cost matrix
    type(cmat) :: c
    
    !dimension of cost matrix
    integer :: dim = 0
    
    !output file unit
    integer :: output = 5
    
    !vector to storing time slots
    integer, allocatable, dimension(:) :: timeslots
    
    !heuristic schedule 
    integer, allocatable, dimension(:) :: heuristicSchedule
    integer, allocatable, dimension(:) :: HSchedule
    
    !statistics related variable (mostly used by branch and bound
    integer:: maxDepth, minDepth, leafCount, validLeafCount, invalidLeafCount, prunedLeafCount 
    real:: avgDepth
    
    !timer related variables
    integer:: tickStart, tickEnd, tickMax, tickRate
    
    contains
    
    !========================================================================
    !initialize ticker 
    !========================================================================
    subroutine initializeTimer()
        implicit none
        call SYSTEM_CLOCK(COUNT_RATE=tickRate, COUNT_MAX=tickMax)
    end subroutine initializeTimer
    
    !========================================================================
    !Start the timer
    !========================================================================
    subroutine startTimer()
        implicit none
        call initializeTimer()
        call SYSTEM_CLOCK(COUNT=tickStart)
    end subroutine startTimer
    
    !========================================================================
    !stop the timer to get the elapsed time
    !========================================================================
    function stopTimer()
        implicit none
        real:: stopTimer, totalTicks
        call SYSTEM_CLOCK(COUNT=tickEnd)
        totalTicks = tickEnd - tickStart
        if(tickEnd < tickStart) then
            totalTicks = totalTicks + tickMax
        end if    
        stopTimer = REAL(totalTicks) / tickRate
    end function stopTimer
    
    !========================================================================
    !prints jobs in the jobtable 
    !========================================================================
    subroutine printjobs()
        implicit none
        integer:: i
        !print the jobtable
        write(output,20) "Index", "|", "p", "|", "d", "|", "w", "|", "r"
        write(output,*) "------------------------------------"
        do i = 1, numjob
            write(output,10) jobtable(i).jobi, "|",  jobtable(i).pi, "|", jobtable(i).di, "|", jobtable(i).wi, "|", jobtable(i).ri
        end do       
        
10      format(I5,2X,A1,2X,I3,2X,A1,2X,I3,2X,A1,2X,I3,2X,A1,2X,I3)
20      format(A5,2X,A1,2X,A3,2X,A1,2X,A3,2X,A1,2X,A3,2X,A1,2X,A3)
    end subroutine printjobs
 
    !========================================================================
    !prints the cost matrix 
    !========================================================================
    subroutine printcostmatrix(nc)
        implicit none
        integer:: i, j
        type(cmat) nc
        write(output,'(A21)') "cost matrix assigned " 
        !print the cost matrix
        write(output,30) 0, (i,i=1,dim)
        write(output,40) "-----", ("--------",i=1,dim)
        do i= 1, dim
            write(output, 30) i, (nc%m(i,j), j=1,dim)   
            !write(*,40) "-----", ("--------",j=1,dim)
        end do
        
30      format(I5,<dim>(1X,"|",1X,I5)) 
40      format(A5,<dim>(A8))    
    
    end subroutine printcostmatrix
    
    !========================================================================
    !prints the cost matrix with row to column assignment marked by *
    !if detail level is 1 then it prints everything
    !if detail level is 2 then it prints only the cost
    !if detail level is 3 then it prints the asssignment 
    !========================================================================
    subroutine printsolvedmatrix(nc, detail)
        implicit none
        integer:: i, j, k
        type(cmat) nc
        integer :: val, detail
        character*1000 row
        character*7 :: cells(dim)
        !print the cost matrix
        write(output, '(A40,I5)') "Cost matrix solved with optimal value = ", nc%z
        if(detail == 3) then
            write(output,50) (i,i=1,dim),(nc%r2c(i),i=1,dim)
        end if    
        if(detail == 1) then 
            write(output,10) 0, (i,i=1,dim)
            write(output,40) "-----", ("--------",i=1,dim)
            do i= 1, dim
                k = 1
                do j =1, dim
                    if (nc%r2c(i)==j) then
                        write(cells(k), 20) nc%m(i,j), '*'
                    else
                        write(cells(k), 20) nc%m(i,j), ' '
                    end if
                    k = k + 1
                end do    
                write(output, 30) i, (cells(j), j=1,dim)   
                !write(*,40) "-----", ("--------",j=1,dim
            end do
        end if
        
10      format(I5,<dim>(1X,"|",1X,I5))       
20      format(1X,I5,A)        
30      format(I5,<dim>("|",A7)) 
40      format(A5,<dim>(A8))   
50      format("mapping =",<dim>(I3,"->"I3,"; "))        
    
    end subroutine printsolvedmatrix
    
    
    !****************************************************************************
    !
    !  SUBROUTINE: readcostmatrix
    !
    !  PURPOSE:  This subroutine reads the cost matrix file and returns the number
    !            of jobs and the cost matrix 
    !  
    !****************************************************************************
    !  INPUT: costfile = cost file path full (see samplecostfile.txt in source folder
    !         for format)        
    !
    !  Output:  costmatrix = cost matrix
    !           n = number of jobs
    !****************************************************************************    
    subroutine readcostmatrix(costfile, n, costmatrix)
        implicit none
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
        implicit none
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
    
    !****************************************************************************
    !
    !  SUBROUTINE: sort job table
    !
    !  PURPOSE:  sorts the jobs by release date and then weight 
    !  
    !****************************************************************************
    !****************************************************************************
    subroutine sortJobs()
        implicit none
        type (jobstruct) :: temp
        INTEGER :: i, j, k, l, ri
        LOGICAL :: swapped = .TRUE.
 
        !sort by bubble sort comparing ready time
        do while (swapped)
            swapped = .FALSE.
            DO i = 1, numjob-1
              !swap if i-th ri is greater then i+1 th ri  
              IF (jobtable(i).ri > jobtable(i+1).ri) THEN
                temp = jobtable(i)
                jobtable(i) = jobtable(i+1)
                jobtable(i+1) = temp
                swapped = .TRUE.    
              END IF
            END DO
        END DO
        
        !now sort by weight
        l = jobtable(1).ri
        j = 1
        do i = 1, numjob
           if(jobtable(i).ri == l) then
               cycle
           else
               if ((i-j) > 1) then
                    !sort by bubble sort comparing weight
                    swapped = .TRUE.
                    DO while (swapped)
                        swapped = .FALSE.
                        DO k = j, i-2
                              !swap if i-th ri is greater then i+1 th ri  
                              IF (jobtable(k).wi > jobtable(k+1).wi) THEN
                                    temp = jobtable(k)
                                    jobtable(k) = jobtable(k+1)
                                    jobtable(k+1) = temp
                                    swapped = .TRUE.    
                              END IF
                        END DO
                    END DO
               end if
                l = jobtable(i).ri
                j = i            
           end if    
        end do    
        
        !reassign job index
        do i=1,numjob
            jobtable(i).jobi = i
        end do    
        
    end subroutine sortJobs
    
    
    !****************************************************************************
    !
    !  SUBROUTINE: assignmatrix
    !
    !  PURPOSE:  this subroutine assign cost matrix c from jobs
    !  
    !****************************************************************************
    !  INPUT: jobfile = full path of job file
    !****************************************************************************
    subroutine assignmatrix()
        implicit none
        integer :: n, i = 1, j = 1, k = 1, p = 1, col = 1, MAX_INTEGER = 32767
        integer :: minRi = 10000 
    
        !calculate the total dimensation fo cost matrix
        ! summation of all processing time as time is granulated as unit time
        do while (i <= numjob)
            dim = dim + jobtable(i).pi
            if (jobtable(i).ri < minRi) then
              minRi = jobtable(i).ri
            end if 
            i = i + 1
        end do    
        
        !allocate the cost matrix
        allocate(c%m(dim,dim))
        allocate(c%r2c(dim))
        allocate(c%c2r(dim))
        allocate(timeslots(dim))
        
        !data cleaning
        !if ri>1 for every i then ri = ri-(rmin-1) because at least one job need to be available in the first slot
        if (minRi > 1) then
            i=1
            do while (i <= numjob)
                jobtable(i).ri = jobtable(i).ri - minRi + 1  
                i = i + 1
            end do
        end if
            
        i = 1
        do while (j <= numjob)
             !every job will contribute its processing time number of rows, each of which denotes one time unit
            do while (p <= jobtable(j).pi) 
                ! populate the colums of row i of cost matrix c
                do while (col <= dim) 
                    !assign a large value to column col of row i of cost matrix c if col < due date of job j or col > dim - processing time of job j
                    if (col < jobtable(j).ri + p -1 .OR. col >= dim - jobtable(j).pi + p + 1)  then
                        c%m(i,col) = MAX_INTEGER                        
                    else
                        !assign 0 last but all parts of job j 
                        if (p < jobtable(j).pi) then
                            c%m(i,col) = 0
                        !for the last part assign the cost of delay calculated by (due date of job j - finishing time) * weight for job j     
                        else 
                            !if (col > jobtable(j).di) then 
                                c%m(i,col) =  jobtable(j).wi * col                                  !(col - jobtable(j).di) * jobtable(j).wi 
                            !else 
                            !    c%m(i,col) = 0
                            !end if
                        end if   
                    end if
                    col = col + 1 !increment col
                end do            
                p = p + 1   !increment p
                i = i + 1   !increment i  
                col = 1 !re-initialize col
                !assign timeslots
                timeslots(k) = j 
                k = k + 1
            end do            
            j = j + 1 !increment j
            p = 1     ! re-initialize p            
        end do      

    end subroutine assignmatrix

    !=====================================================================
    ! returns false is the assignment is valid i.e no blocked cell is selected
    !=====================================================================
    function isAssignmentVaild(sol)
        type(solution) sol
        logical :: isAssignmentVaild, isValid
        
        integer, dimension(dim) :: r2c
        integer, dimension(dim,dim) :: cost
        integer :: i
        
        cost(1:dim,1:dim) = sol%c%m(1:dim,1:dim)
        r2c(1:dim) = sol%c%r2c(1:dim)
        isValid = .TRUE.
        !check if any of the blocked cell is selected 
        do i = 1 , dim
            if (cost(i,r2c(i)) > 9999) then
                isValid = .FALSE.
            end if
        end do    
        isAssignmentVaild = isValid
    end function
    
    !======================================================================================
    !Implementation of WSRPT heuristics 
    !Batsyn, Goldengorin, Pardalos, Sukhov 2013
    !======================================================================================
    subroutine wsrpt(optimumCost, schedule)
        
        implicit none
        
        integer:: optimumCost
        integer, dimension(dim):: schedule
        integer, dimension(numJob):: rJobs
        integer:: i,j,k,l
        real :: bestRatio = 0.0, ratio, w, p
        logical:: canAssign = .FALSE.
        
        !populate the remaining processing time to rJobs
        do j=1,numJob
            rJobs(j) = jobtable(j).pi
        end do    
        
        !initialize schedule
        do j=1,dim
            schedule(j) = 0
        end do
        
        do i=1,dim
            !collect all assignable jobs 
            canAssign = .FALSE.
            do j = 1,numjob
                if(jobtable(j).ri <= i .AND. rJobs(j) > 0) then
                     w = jobtable(j).wi 
                     p = rjobs(j)
                     ratio = w/p
                     if(ratio > bestRatio) then
                         canAssign = .TRUE.
                         bestRatio = ratio
                         k = j
                     end if    
                end if   
            end do
            !assign the job with best wi/rhoi ratio and decrease the remaining processing time
            if (canAssign) then
                schedule(i) = k
                rJobs(k) = rJobs(k) - 1
            end if
            bestRatio = 0.0
        end do 
        
        do j=1,numJob
            rJobs(j) = jobtable(j).pi
        end do 
        
        !calculate optimum cost 
        optimumCost = 0
        do i=1,dim
            if (schedule(i) > 0) then 
                rJobs(schedule(i)) = rJobs(schedule(i)) - 1
                if (rJobs(schedule(i)) == 0) then
                    optimumCost = optimumCost + jobtable(schedule(i)).wi * i
                end if  
            end if  
        end do
        
        !asign HSchedule (schedule by part to match matrix formation technique)
        allocate(HSchedule(dim))
        l = 1
        iloop: do i=1,numjob
            kloop: do k = 1,jobtable(i).pi
                jloop: do j=1,dim
                    if (schedule(j)==i) then
                        schedule(j) = 0
                        HSchedule(l) = j
                        l = l + 1
                        EXIT jloop
                    end if   
                end do jloop   
            end do kloop
        end do iloop   
        
        write(output, '(A40,I5)') "Optimum value of the schedule by WSRPT = ", optimumCost
        write(output, 10) (HSchedule(i),i=1,dim) 
        
10      format("{",<dim>(I3," "),"}")        
    
    end subroutine
    
end module    