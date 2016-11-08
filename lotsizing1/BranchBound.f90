module BranchBound
    
    !type matpos
    !    integer :: i
    !    integer :: j
    !end type matpos
    
    type violation
        integer ct(2)
        integer nct(2)
    end type violation
    
    type(violation) vios(100)
    
    integer :: violationSize = 0, upperBound = 0
    
    contains   
    
    !-----------------------------------------------------------
    !main subroutine for branch and bound
    !-----------------------------------------------------------
    subroutine DFSBB(RC)
        use StackObject
        use global
        use lapjv
        implicit none
        
        type(cmat) RC
        
        integer :: x(dim), &                !col assigned to row
                   y(dim), &                !row assigned to column
                   u(dim), &                !Dual Row variable
                   v(dim), &                !Dual column variable
                   z, &    
                   KK(2500), &
                   FIRST(101)
        
        !declare stack to store nodes to be visited
        !size of stack is fixed to 10, can be increased later
        type(StackT) s 
        s = NewStack(10)
        
        !initialize the stack with the root element
        call push(RC, s)
        
        !continue looping until there is no more nodes to visit
        do while(.not. isStackEmpty(s))
           
           !get the first node from to be visited list
           !the stack will always return a sub problem
           RC = top(s)
           call pop(s)
           
           !run the cost matrix through LAPJV
           call JOVOFD(dim, RC%m, x, y, u, v, z)
           !call JOVOSAP(dim,RC%m,KK,FIRST,X,Y,U,V,z)
           call printsolvedmatrix(RC, y, x, z)           
           
           !identify every violations and classify them in completion and non completion time groups
           !call classifyViolations(RC, x, y) 
           
           !enumerate each violation and calculate bottleneck upper tolerance 
           
            
        enddo
        
    
    end subroutine DFSBB
    
    !this subroutine classify violations into completion and non-completions pairs 
    !this subroutine populates the vios array with completion and non-compeltion pairs
    subroutine classifyViolations(nc, x, y)
        use global    
        implicit none
    
        type(cmat) nc
        integer:: x(dim), & !column assigned row 
                  y(dim)    !row assigned to column
        
        integer :: i = 1, j = 1, k = 1, firstRow = 0, lastRow = 0
        type(violation) v
        
        !for every job
        do i = 1, numjob
            !calculate the first and last row (timeslot) dedicated to this job
            firstRow = lastRow + 1
            lastRow = lastRow + jobTable(i).pi            
            !loop over only non-completion rows
            do j = firstRow, lastRow-1
                !if col assigned to non-completion row is greater than the col assigned to compeltion time row
                if (x(lastRow) < x(j)) then
                    !mark as violation
                    v.ct(1) = lastRow
                    v.ct(2) = x(lastRow)
                    v.nct(1) = j
                    v.nct(2) = x(j)
                    vios(k) = v
                    k = k + 1
                end if
            end do    
        end do    
        violationSize = k
        
    end subroutine classifyViolations
    
    !this subroutine updates the bottleneck tolerance 
    subroutine calculateBottleneckTolerance(nc, upperBound)
        use global
        implicit none
        integer :: upperBound, i = 1, j = 1 , k = 1, minCT = 0, minNCT = 0
        type(cmat) nc
        !loop though all violations
        do i = 1,  violationSize
           !calculate change if ct violation is excluded
           j =  solveJVByExcluding(nc, vios(i).ct)
           if (minCT < (j - upperBound)) then
                minCT = j - upperBound
           end if
           !calculate change if nct violation is excluded
           j =  solveJVByExcluding(nc, vios(i).nct)
           if (minNCT < (j - upperBound)) then
                minNCT = j - upperBound
           end if
        end do 
        return max(minCT, minNCT)
    
    end subroutine calculateBottleneckTolerance
    
    function solveJVByExcluding(nc, matpos)
        use global
        implicit none
        integer :: matpos(2), solveJVByExcluding        
        type(cmat) nc
        integer :: large, &
                   x(dim), &                !col assigned to row
                   y(dim), &                !row assigned to column
                   u(dim), &                !Dual Row variable
                   v(dim), &                !Dual column variable
                   z         
        
        !block the cell
        nc.m(matpos(1),matpos(2)) = huge(large)
        
        call JOVOFDTEST(dim, nc.m, x, y, u, v, z)
        
        solveJVByExcluding = z        
    
    end function solveJVByExcluding
    
    
end module BranchBound    