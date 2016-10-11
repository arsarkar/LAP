module BranchBound
    
    type matpos
        integer :: i
        integer :: j
    end type matpos
    
    type violation
        type(matpos) ct
        type(matpos) nct
    end type violation
    
    type(violation) vios(50)
    
    contains
    
    
    
    !-----------------------------------------------------------
    !main subroutine for branch and bound
    !-----------------------------------------------------------
    subroutine DFSBB(RC)
        use StackObject
        use global
        implicit none
        
        type(cmat) RC
        
        integer :: x(dim), &                !col assigned to row
                   y(dim), &                !row assigned to column
                   u(dim), &                !Dual Row variable
                   v(dim), &                !Dual column variable
                   z         
        
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
           call JOVOFDTEST(dim, RC%m, x, y, u, v, z)
           call printsolvedmatrix(RC, x, y)           
           
           !identify every violations and classify them in completion and non completion time groups
           
            
        enddo
        
    
    end subroutine DFSBB
    
    subroutine classifyViolations(nc)
        use global    
        implicit none
    
        type(cmat) nc
        
        integer :: i = 1
        
        do while (i <= numJob)
            
            
            
        end do    
        
    end subroutine classifyViolations
    
    
    
end module BranchBound    