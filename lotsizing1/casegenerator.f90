module casegenerator
    use global
    
    !job table
    type (jobstruct), dimension(:), allocatable :: jobs
    logical :: randSeeded = .FALSE.
    integer :: MAX_PI = 100, MAX_WI = 20, MIN_WI = 1
    integer, dimension(:), allocatable :: RANDSEED
    double precision :: relRangeDueDate, avgTardFact
    
    contains
    
    !========================================================================
    ! call this subroutine to change the following default values
    ! MAX_PI = 100, MAX_WI = 20, MIN_WI = 1
    !========================================================================
    subroutine configure(maxpi, maxwi, minwi, seed)
        implicit none    
        
        integer :: maxpi, maxwi, minwi, seed
        
        MAX_PI = maxpi
        MAX_WI = maxwi
        MIN_WI = minwi
        allocate(RANDSEED(1))
        RANDSEED(1) = seed
        
    end subroutine configure
    
    !========================================================================
    ! Random number generator; generates integer random number between min and max    
    !========================================================================
    function rand(min, max)
        integer :: rand, min, max
        double precision :: r
                
        !seed the random number generator
        if(.not. randSeeded) then
            call RANDOM_SEED(PUT=RANDSEED)       
            randSeeded = .TRUE.
        end if
        
        !generate random number to scale it to the max and min
        call RANDOM_NUMBER(r)
        
        rand = (r * (max + 1 - min)) + min
        
    end function rand

    !========================================================================
    !this subroutine generates jobs and fills up jobs table 
    ! this subroutine follows the strategy devised by Potts and Wassenhove (1982)
    ! and reimplemented by Bulbul et al (2007)
    ! input::
    ! n = number of jobs 
    ! r = relative range of due dates U[0.2,1.0]
    ! t  = Avg. tardiness factor U[0.2,0.8]
    !========================================================================    
    subroutine potts1982(n, r, t)
        use global
        implicit none
        
        integer:: n, i, P = 0, C = 0
        double precision :: r, t
        type (jobstruct) :: job
        
        !redim the array of jobs to number of jobs
        allocate(jobs(n))
        relRangeDueDate = r
        avgTardFact = t
        
        !assign Pi from U[1,100]
        do i = 1, n
           job%jobi = i
           job%pi = rand(1,MAX_PI)
           P = P + job%pi
           jobs(i) = job
        end do    
    
        !assign Di, Ri and Wi
        do i = 1, n
           C  = C + jobs(i)%pi
           jobs(i)%di = rand(ceiling(P*(1-t-r/2)), ceiling(P*(1-t+r/2)))
           jobs(i)%ri = rand(0,C/2)
           jobs(i)%wi = rand(MIN_WI, MAX_WI)
        end do  
        
    end subroutine potts1982
    
    
    !function writeJobs(filePath)
    !    
    !    !character(150) :: filePath
    !    !integer :: jobfileUnit = 20, iostat
    !    !
    !    !!!open the file to write
    !    !open(unit = jobfileUnit, file = filePath, iostat=iostat)
    !    !
    !    !if (iostat == 0) then
    !    !!!write the file 
    !    !write(jobfileUnit, 10 
    !    !
    !    !end if     
    !    
    !end function
    

end module casegenerator    