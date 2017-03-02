module single_premp_p_r_wt
use global
use casegenerator
use MHEAP
implicit none
    
    type(THEAP) :: heap
    type(jobstruct), dimension(:), allocatable:: sJobs
    integer, dimension(:), allocatable:: schedule
    
    contains   

    !Subroutine SlackX1 heuristic 
    subroutine slackX1(tt)
        use global
        implicit none 
        integer :: tt, i = 1, fJobs(numjob), t = 0
        double precision:: node(1)
        
        !initialize the job and schedule array
        allocate(sJobs(numjob))
        sjobs(1:numjob) = jobtable(1:numjob)
        
        !calculate total length of horizon
        do while (i <= numjob)
            dim = dim + sJobs(i).pi
            i = i + 1
        end do        
        allocate(schedule(dim))
        
        !initialize the heap
        call heap%INIT(numjob, 1, job_cmp)
    
        !populate the heap with jobs
        do i = 1, numjob
            call heap%INSERT([DBLE(i)])
        end do
        
        !populate the schedule
        do i = 1, dim
            !pop the root node, max heap should return the job with maximum slackness 
            call heap%POP(node)
            !schedule the job 
            schedule(i) = int(node(1))
            !decrease the processing time by one unit
            sJobs(schedule(i)).pi = sJobs(schedule(i)).pi - 1
            !insert the job back again if it still remaining to be processes
            if (sJobs(schedule(i)).pi > 0) then
                call heap%INSERT(node)
            else
                !mark the finishing time
                fJobs(schedule(i)) = i
            end if
        end do    
                
        !calculate total tardiness
        do i = 1, numjob
            !get the finishing time
            t = t + max(0, sjobs(i).di - fJobs(i))
        end do    
        tt = t
        
    end subroutine slackX1

    !comparison function for max heap
    logical function job_cmp(j1, j2)
        double precision, intent(in) :: j1(:), j2(:)
        job_cmp = getSlack(j1(1)) > getSlack(j2(1))        
    end function job_cmp
    
    !calculate the slack for job index j
    double precision function getSlack(j)
    implicit none
        double precision:: j
        getSlack = max(0, (sJobs(int(j)).di - (sJobs(int(j)).ri - sJobs(int(j)).pi)))
    end function getSlack
    
    
end module     