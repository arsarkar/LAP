module single_premp_p_r_wt
use global
use casegenerator
use MHEAP
implicit none
    
    type(THEAP) :: heap
    type(jobstruct), dimension(:), allocatable:: sJobs

    contains
    
    

    subroutine slackX1(schedule)
        use global
        implicit none 

        integer :: horizonLength, i
        integer, dimension(dim):: schedule
        
        horizonLength = 0
        !do i = 1 , dim
        !    horizonLengh = horizonLengh + jobTable(i).pi
        !end do    
        
        
        !call heap%INIT(2, 2, job_cmp)
    
    
    end subroutine slackX1


    logical function job_cmp(j1, j2)
        double precision, intent(in) :: j1, j2
        job_cmp = getSlack(j1) < getSlack(j2)        
    end function job_cmp
    
    !calculate the slack for job index j
    double precision function getSlack(j)
    implicit none
        double precision:: j
        
        
        getSlack = 0
    end function getSlack
    
end module     