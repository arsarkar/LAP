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
    
    contains
    
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
    subroutine printsolvedmatrix(nc, x, y, z, detail)
        implicit none
        integer:: i, j, k, z
        type(cmat) nc
        integer :: x(dim), y(dim), val, detail
        character*1000 row
        character*7 :: cells(dim)
        !print the cost matrix
        write(output, '(A40,I5)') "Cost matrix solved with optimal value = ", z
        if(detail == 3) then
            write(output,50) (i,i=1,dim),(x(i),i=1,dim)
        end if    
        if(detail == 1) then 
            write(output,10) 0, (i,i=1,dim)
            write(output,40) "-----", ("--------",i=1,dim)
            do i= 1, dim
                k = 1
                do j =1, dim
                    if (x(i)==j) then
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
    
        !calculate the total dimensation fo cost matrix
        ! summation of all processing time as time is granulated as unit time
        do while (i <= numjob)
            dim = dim + jobtable(i).pi
            i = i + 1
        end do    
        
        !allocate the cost matrix
        allocate(c%m(dim,dim))
        allocate(c%r2c(dim))
        allocate(c%c2r(dim))
        allocate(timeslots(dim))

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
        integer:: i,j,k
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
        do i=1,dim
            if (schedule(i) > 0) then 
                rJobs(schedule(i)) = rJobs(schedule(i)) - 1
                if (rJobs(schedule(i)) == 0) then
                    optimumCost = optimumCost + jobtable(schedule(i)).wi * i
                end if  
            end if  
        end do 
        
        write(output, '(A40,I5)') "Optimum value of the schedule by WSRPT = ", optimumCost
        write(output, 10) (schedule(i),i=1,dim) 
        
10      format("{",<dim>(I3," "),"}")        
    
    end subroutine
    
end module    