module branchbound
    use global
    use StackObject
    
    type(violation) vios(100)
    
    integer :: violationSize = 0, upperBound = 9999, reportInterval = 100

    !declare stack to store nodes to be visited
    !size of stack is fixed to 10, can be increased later
    type(StackT) unvisited     
    
    contains   
    
    !-----------------------------------------------------------
    !main subroutine for branch and bound
    !-----------------------------------------------------------
    subroutine DFSBB(RC)
        use StackObject
        use global
        use jv
        implicit none
        
        type(cmat) RC
        
        integer :: costmat(100,100), &
                   x(dim), &                !col assigned to row
                   y(dim), &                !row assigned to column
                   u(dim), &                !Dual Row variable
                   v(dim), &                !Dual column variable
                   z, &    
                   KK(2500), &
                   FIRST(101), &
                   large, &
                   solutionCount, &
                   tempCount, &
                   stackSize
        Logical isValid, hasLowerBound, hasUpperBound
        REAL :: elapsedTime, timeRemaining
        
        !placeholder for solution
        type(solution) sol, soltmp
        type(solution) solx,soln, solfinal
        unvisited = NewStack(100)
        
        !initialize the stack with the root element and other counters/flags 
        sol.c = RC
        solutionCount = 0
        sol.depth = 0
        maxDepth = 0
        minDepth = 9999
        avgDepth = 0
        leafCount = 0
        validLeafCount = 0
        invalidLeafCount = 0
        prunedLeafCount = 0
        tempCount = 0
        
        !write the heuristic solution assigned to the cost matrix
        sol%c%r2c = Hschedule
        sol%c%z = upperBound
        call printsolvedmatrix(sol.c, 1)
        
        !solve the solution for the first time so that 
        !the stack always stores solved solutions
        !costmat(1:dim,1:dim) = sol.c.m(1:dim,1:dim)
        !run the cost matrix through LAPJV
        !call JOVOFD(dim, costmat, x, y, u, v, z)
        call solveLAP(sol.c.m, x, y, z)
        sol%c%r2c = x
        sol%c%c2r = y
        sol%c%z = z

        !call printsolvedmatrix(sol.c, 1) 
        !call printsolvedmatrix(sol.c, 3) 
        
        call push(sol, unvisited)
        
        write(output, 50)
        
        ! start the timer
        call startTimer()
        
        !continue looping until there is no more nodes to visit
        do while(.not. isStackEmpty(unvisited))
           
           write(output, 50)
           write(output, 10) 
           write(output, 50)
           
           !get the first node from to be visited list
           !the stack will always return a sub problem
           sol = top(unvisited)
           call pop(unvisited)
           soltmp = sol 
           solutionCount = solutionCount + 1 
           tempCount = tempCount + 1
           
           write(output, 30) "SolutionCount :: ", solutionCount
           write(output, 30) "Depth of solution :: ", sol%depth
           if(sol.sr == 1) then
               write(output, 70) "including ", sol%pivot(1), sol%pivot(2) 
           else if(sol.sr == 2) then
               write(output, 70) "excluding ", sol%pivot(1), sol%pivot(2) 
           end if    
           write(output, 30) "Current upper bound =", upperBound
           call printsolvedmatrix(sol.c, 3)
           
           if (tempCount > reportInterval) then
               elapsedTime = stopTimer()
               stackSize = getCount(unvisited)
               timeRemaining = (stackSize * (elapsedTime / solutionCount)) / 60
               tempCount = 0
               write(*, '(A50)') "----------------------------------------------------------------"
               write(*, 30) "SolutionCount :: ", solutionCount
               write(*, 30) "Depth of solution :: ", sol%depth
               write(*, 30) "Current upper bound =", upperBound
               write(*, 30) "Solution cost =", sol%c%z
               write(*, 30) "Estimated time remaining (min) :: ", stackSize
           end if
           
           !calculate depth statistics
           if(sol%depth > maxDepth) then
               maxDepth = sol%depth
           end if
           
           !prune the tree by checking the upperBound bound
           if (sol%c%z > upperBound) then
               write(output, 30) "Pruning the tree with solution = ", sol%c%z
               leafCount = leafCount + 1
               prunedLeafCount = prunedLeafCount + 1
               avgDepth = avgDepth + sol%depth
               if(sol%depth < minDepth) then
                   minDepth = sol%depth
               end if
               cycle
           end if    
           
           !check if the solution is valid, i.e no blocked cell is selected (this can be checked using upper bound, when blocked cell is selected it will surely cross upper bound)
           isValid = isAssignmentVaild(sol)
           if (.NOT. isValid) then
               write(output, 30) "Pruning the tree with invalid solution at ", solutionCount
               leafCount = leafCount + 1
               invalidLeafCount = invalidLeafCount  + 1
               avgDepth = avgDepth + sol%depth
               if(sol%depth < minDepth) then
                   minDepth = sol%depth
               end if               
               cycle
           end if  
           
           write(output, 50)
           
           !check if the solution has sequence constrin maintained, if not then populate every violation in vios array by
           !identifying every violations and classifying them in completion and non completion time groups
           !if a valid solution is found then cost of the solution is recorded as lower bound and 
           !go to the next solution
           isValid = classifyViolations(sol)
           if(isValid) then
                !for a valid solution lower bound is recorded, Can a solution go lower than this? probably not being a feasible at the same time
                if (sol%c%z <= upperBound) then
                    upperBound = sol%c%z
                    solfinal = sol
                end if
                write(output, 30) "Solution is valid! new upperBound = ", upperBound
                leafCount = leafCount + 1
                validLeafCount = validLeafCount + 1
                avgDepth = avgDepth + sol%depth
               if(sol%depth < minDepth) then
                   minDepth = sol%depth
               end if
               cycle
           end if    
           
           write(output, 50)
           
           !enumerate each violation and calculate bottleneck upper tolerance 
           call calculateBottleneckTolerance(sol, 1)         
           
           write(output, 50)
           
           !the maximum solution found is the sub problem which is already solved by ~excluding~ the pivot
           solx = soltmp
           !here we generate the solution by including the pivot
           soln = soltmp
           call generateSubProblem(solx,soln)
        
        end do
        
        !stop the time
        elapsedTime = stopTImer()
        avgDepth = avgDepth / leafCount
        
        write(output,50)
        write(output,50)
        write(output,50)
        
        write(output, 40) solfinal%c%z
        call printsolvedmatrix(solfinal.c, 1)
        write(output,50)
        write(output,50)
        write(output, 30) "Total nodes = ", solutionCount
        write(output, 30) "Maximum depth of solution tree = ", maxDepth
        write(output, 30) "Minimum depth of solution tree = ", minDepth
        write(output, 60) "Average depth of solution tree = ", avgDepth
        write(output, 30) "Total number of leaves in solution tree = ", leafCount
        write(output, 30) "Total number of valid leaves in solution tree = ", validLeafCount
        write(output, 30) "Total number of invalid leaves in solution tree = ", invalidLeafCount
        write(output, 30) "Total number of pruned leaves in solution tree = ", prunedLeafCount
        write(output, 60) "Total time taken = ", elapsedTime

10      format("----------------------solving new subproblem-------------------------")
20      format("Solving subproblem by excluding (",I3,",",I3,") with upper bound ")
30      format(A50, I10)        
40      format("Optimum cost = ", I3, " optimum assignment ->")       
50      format("------------------------------------------------------------------------------------------")   
60      format(A50, F5.2)    
70      format(A50, "(",I3,",",I3,")")        
    
    end subroutine DFSBB
    
    !========================================================================================
    !this subroutine classify violations into completion and non-completions pairs 
    !this subroutine populates the vios array with completion and non-compeltion pairs
    !========================================================================================
    function classifyViolations(sol)
        use global    
        implicit none
    
        type(solution) sol
        integer:: x(dim), & !column assigned row 
                  y(dim)    !row assigned to column
        logical:: classifyViolations
        
        integer :: i = 1, j = 1, k = 0, firstRow = 0, lastRow = 0
        type(violation) v 
        
        x = sol%c%r2c
        y = sol%c%c2r
        lastRow = 0
        k = 0
        !write(output,30)
        !call printsolvedmatrix(sol%c, sol%c%r2c, sol%c%c2r, sol%C%z, 1)
        
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
                    k = k + 1
                    !globally add the violation to violation collection
                    vios(k) = v
                end if
            end do    
        end do    
        !globally assign violation size 
        violationSize = k
        
        !print all violation
        write(output,20) violationSize
        !write(output,10) (vios(i).nct(1),vios(i).nct(2),vios(i).ct(1),vios(i).ct(2),i=1,violationSize) 
        
        if(violationSize .EQ. 0) then
            classifyViolations = .TRUE.
        else
            classifyViolations = .FALSE.
        endif    
        
10      format(<violationSize>("(",I3,",",I3,")-(",I3,",",I3,");")) 
20      format("total number violations found",I3)        
30      format("classifying violations for cost matrix....") 
        
    end function classifyViolations
    
    !===================================================================================
    !this subroutine updates the bottleneck tolerance 
    !===================================================================================
    ! sol ->  solution with cost matrix but no pivot and branching rule yet decided
    ! strategy -> different ways of calculating bottleneck tolerance 1->min-min-max 2->max-ct
    !===================================================================================
    subroutine calculateBottleneckTolerance(sol, strategy)
        use global
        use StackObject
        implicit none
        integer :: upperBound, i = 1, z = 1 , k = 1, l, minCT = 0, minNCT = 0, maxCT = 0, bottleNeck, strategy
        integer, dimension(2):: minCTIndex, maxCTIndex, minNCTIndex, maxProbIndex 
        type(solution) sol
        type(solution) sol0
        type(solution) minCTSol, maxCTSol, minNCTSol
        type(solution) subprobs(100)
        
        sol0 = sol
        minCT = 0
        minNCT = 0
        k=1
        !loop though all violations
        do i = 1,  violationSize
           
           !calculate change if ct violation is excluded
           !write(output,20) vios(i).ct(1), vios(i).ct(2)
           sol = sol0                                                   !get the parent solution
           !generate new problem pivoting on ct vviolation
           call calculateTolerance(sol, vios(i).ct)
           subprobs(k) = sol
           k = k + 1
            
           !calculate the minimum ct tolerance 
           if (i==1) then
               minCT = sol.tolerance
               minCTSol = sol
               minCTIndex = vios(i).ct
               maxCT = sol.tolerance
               maxCTSol = sol
               maxCTIndex = vios(i).ct
           else    
               if (sol.tolerance < minCT) then
                    minCT = sol.tolerance
                    minCTSol = sol
                    minCTIndex = vios(i).ct
               end if
               if (sol.tolerance > maxCT) then
                   maxCT = sol.tolerance
                   maxCTSol = sol
                   maxCTIndex = vios(i).ct
               end if
           end if
           
           !!calculate change if nct violation is excluded           
           !write(output,10) vios(i).nct(1), vios(i).nct(2)
           !sol = sol0  
           !
           !!generate new problem based on nct violation
           !call calculateTolerance(sol, vios(i).nct)
           !!put the new solution in the ascending order of tolerance
           !subprobs(k) =  sol 
           !k = k +1 
           ! 
           !!calculate the minimum nct tolerance
           !if (i==1) then
           !    minNCT = sol.tolerance
           !    minNCTSol = sol
           !    minNCTIndex = vios(i).nct
           !else   
           !    if (sol.tolerance < minNCT) then
           !         minNCT = sol.tolerance
           !         minNCTSol = sol
           !         minNCTIndex = vios(i).nct
           !    end if
           !end if
            
        end do 
        
        !calculate bottle nexk tolerance 
        !set the sol with the min-min-max solution
        !set the pivot 
        !if (minCTSol.tolerance >= minNCTSol.tolerance) then
            bottleNeck = maxCT
            sol = maxCTSol
            sol%pivot = maxCTIndex
        !else
        !    bottleNeck = minNCT
        !    sol = minNCTSol
        !    sol%pivot = minNCTIndex
        !end if
        
        write(output,30) sol%pivot(1), sol%pivot(2), bottleNeck
        
10      format("Solving with excluded NCT violation (",I3,",",I3,")")  
20      format("Solving with excluded CT violation (",I3,",",I3,")")  
30      format("Branching pivot found (",I3,",",I3,") with bottleneck tolerance ",I5)         
        
    end subroutine calculateBottleneckTolerance
    
    !===========================================================
    !generate a subproblem by pivoting on the violation     !
    !===========================================================
    subroutine calculateTolerance(sol, violation)
        use global
        implicit none
        type(solution) sol
        integer :: violation(2), k,  i, z, t, LB
        
        !lower bound is the solution 
        LB = sol%c%z
        !calculate change if ct violation is excluded                                                  
        !get the parent solution
        sol.pivot = violation                                       
        !set the pivot
        !solve the new solution by excluding the pivot
        z =  solveJVByExcluding(sol, violation, 1) 
        t = z - LB                                           !calculate tolerance
        sol.tolerance = t                                    !set the tolerance to solution
        
    end subroutine calculateTolerance
    
    !=============================================================
    !Solve the cost matrix by blocking the matpos supplied
    !=============================================================
    function solveJVByExcluding(nc, matpos, strategy)
        use global
        use jv
        implicit none
        integer :: matpos(2), solveJVByExcluding, strategy, i, j, k, l       
        type(solution) nc
        integer :: COSTMAT(100,100), &
                   large = 10000, &
                   x(dim), &                !col assigned to row
                   y(dim), &                !row assigned to column
                   u(dim), &                !Dual Row variable
                   v(dim), &                !Dual column variable
                   z
        
        write(output, 10) matpos(1), matpos(2), strategy
        !block the cell
        if (strategy==1) then 
            nc.c.m(matpos(1),matpos(2)) = large
        else if (strategy==2) then
            j = timeslots(matpos(1))
            k = matpos(2) - jobtable(j).pi + 1
            do i = matpos(1)-jobtable(j).pi+1 , matpos(1)
                do l = k , dim
                    nc.c.m(i,l) = large
                end do
                k = k + 1
            end do    
        end if
        
        call solveLAP(nc.c.m, x, y, z)
        NC%sr = 2
        NC%c%z = z
        NC%c%r2c = x
        NC%c%c2r = y
        call printsolvedmatrix(NC.c, 1) 
        
        solveJVByExcluding = z 
        
10      format("Solving cost matrix y excluding (",I3,",",I3,") with strategy-",I3 )        
    
    end function solveJVByExcluding    
        
    !=======================================================================================
    !takes solution by exclusding pivot as input and produces solution by including pivot
    !both of these solutions is added to the stack
    !=======================================================================================
    subroutine generateSubProblem(solx,soln)
        
        type(solution) solx,soln
        integer :: cost
        
        !change the pivots to assignments as per the heuristic schedule
        solx%pivot(2) = HSchedule(solx%pivot(1))
        
        cost = solveJVByIncluding(soln, solx%pivot, 4)
        cost = solveJVByExcluding(solx, solx%pivot, 2)
        
        !increase the depth of solutions at branch
        solx%depth = solx%depth + 1
        soln%depth = soln%depth + 1
        
        !should we first choose including or excluding? 
        !here whatever is minimum is chosen first
        if(solx.c.z >= soln.c.z) then
            call push(solx, unvisited)  !will chosen later
            call push(soln, unvisited)
        else
            call push(soln, unvisited)  !will chosen later
            call push(solx, unvisited)
        end if
    
    end subroutine generateSubProblem
    
    !=======================================================================================
    !Solve the cost matrix by including the matpos supplied 
    !this is done by two strategies 
    !strategy 1->horizontal blocking 2->vertical blocking 3-> mixed blocking
    !=======================================================================================
    function solveJVByIncluding(nc, matpos, strategy)
        use global
        use jv
        implicit none
        integer :: matpos(2), solveJVByIncluding, strategy        
        type(solution) nc
        integer :: COSTMAT(100,100), &
                   nccopy(dim,dim), & 
                   large = 10000, &
                   x(dim), &                !col assigned to row
                   y(dim), &                !row assigned to column
                   u(dim), &                !Dual Row variable
                   v(dim), &                !Dual column variable
                   ass(dim), &
                   z, &
                   i, &
                   j, &
                   k, &
                   l
        
        write(output, 10) matpos(1), matpos(2), strategy
        !call printcostmatrix(NC%c)
        nccopy(1:dim,1:dim) = NC%c%m(1:dim,1:dim)
        ass(1:dim) = timeslots(1:dim)
        !block the cells to include the cell matpos
        if (strategy==1) then
            do i=1, dim
                !block the CT row except the cell to be included 
                if(i /= matpos(2)) then
                    nccopy(matpos(1),i) = large
                end if
                !block the NCT rows
                j = ass(matpos(1))
                k = matpos(1) - 1
                l = matpos(2)
                do while(j == ass(k))
                    if(i >= l) then 
                        nccopy(k,i) = large
                    end if
                    k = k - 1
                    if(k == 0) then
                        exit
                    end if    
                    l = l - 1
                end do    
            end do 
        else if(strategy==2) then
            do i=1, dim
                if(i /= matpos(1)) then
                    nccopy(i,matpos(2)) = large
                end if
            end do             
        else if(strategy==3) then
            do i=1, dim
                if(i /= matpos(2)) then
                    nccopy(matpos(1),i) = large
                end if
            end do
            do i=1, dim
                if(i /= matpos(1)) then
                    nccopy(i,matpos(2)) = large
                end if
            end do 
        else if(strategy==4) then
            do i=1 , matpos(2)-1
                nccopy(matpos(1),i) = large
            end do    
            j = ass(matpos(1))
            i = matpos(1)-1
            l=matpos(2)+1
            do while(j == ass(i))
                do k = l , dim
                    nccopy(i,k) = large
                end do
                i = i - 1
            end do
        end if    
        
        !costmat(1:dim,1:dim) = nccopy(1:dim,1:dim)
        !call printcostmatrix(NC%c)
        !call JOVOFD(dim, costmat, x, y, u, v, z)        
        call solveLAP(nccopy, x, y, z)
        nc%c%m(1:dim,1:dim) = nccopy(1:dim,1:dim)
        nc%sr = 1
        nc%c%z = z
        nc%c%r2c = x
        nc%c%c2r = y  
        nc%pivot = matpos
        call printsolvedmatrix(NC.c, 1)  
        
        solveJVByIncluding = z      

10      format("Solving cost matrix y including (",I3,",",I3,") with strategy ",I3)         
    
    end function solveJVByIncluding
    
end module BranchBound    