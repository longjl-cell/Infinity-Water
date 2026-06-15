
   module GAUSS_BASIC_MODULE
    
    implicit none
    
    type GAUSS_BASIC_TYPE    
        
        integer, PUBLIC                              ::  NUMBER_POINTS  = 10
        real(8),  PUBLIC, dimension(:), ALLOCATABLE  ::  SEEDS_POINT 
        real(8),  PUBLIC,  dimension(:), ALLOCATABLE ::  WEIGHTS_POINT 
        LOGICAL,  PUBLIC                             ::  SEEDS_FORMING_FLAG  = .FALSE.
        
    contains
   
    procedure   ::   ALLCOATE_ARRAY
    procedure   ::   CALCUATE_SEEDS_WEIGHTS
    
    
    end type GAUSS_BASIC_TYPE 
    
    contains
    
    subroutine ALLCOATE_ARRAY( THIS, NUMBER_POINTS )
    class(GAUSS_BASIC_TYPE), intent( inout) :: THIS
    integer,   optional,   intent( in ) ::  NUMBER_POINTS
    
 !-----------------LOCAL VARIABLE---------------------!   
    
    integer  ::    ALLCOATE_ARRAY_FLAG
    
    
   if( PRESENT(NUMBER_POINTS )  )then
       THIS%NUMBER_POINTS = NUMBER_POINTS 
   end if
   
   
    
    if( THIS%SEEDS_FORMING_FLAG == .TRUE.   )then
  
       deallocate(THIS%SEEDS_POINT,  THIS%WEIGHTS_POINT,   stat =  ALLCOATE_ARRAY_FLAG )
    
           if( ALLCOATE_ARRAY_FLAG /=0 )then
               
              write(*,*)  'ERROR IN DEALLOCATE ERROR GAUSS INTEGRATION METHOD' 
              PAUSE
              
           end if 
 
    end if
    
    allocate( THIS%SEEDS_POINT( THIS%NUMBER_POINTS),   THIS%WEIGHTS_POINT( THIS%NUMBER_POINTS), &
            stat =  ALLCOATE_ARRAY_FLAG )
    
    if( ALLCOATE_ARRAY_FLAG /=0 )then
        
       write(*,*)  'ALLOCATE ERROR GAUSS INTEGRATION METHOD' 
       PAUSE
       
    end if
    
    THIS%SEEDS_POINT(:) = 0.0D0
    THIS%WEIGHTS_POINT(:) = 0.0D0
    
    THIS%SEEDS_FORMING_FLAG = .TRUE. 
    
    
    end  subroutine ALLCOATE_ARRAY
    
    
    
    subroutine CALCUATE_SEEDS_WEIGHTS( THIS  )
    class(GAUSS_BASIC_TYPE), intent( inout) :: THIS
  !  integer,      intent( in ) ::  NUMBER_POINTS
    
    
    end  subroutine CALCUATE_SEEDS_WEIGHTS
   
    
 
    
    
    
    
 end module GAUSS_BASIC_MODULE    
    
    
    
    
    
 module GAUSS_LEDENDRE_MODULE
    
    use GAUSS_BASIC_MODULE
    use precision_mod
 
    IMPLICIT NONE
    
    type, extends(GAUSS_BASIC_TYPE) :: GAUSS_LEGENDRE_TYPE
   
        integer   ::   NEWTON_ITEARTION_MAX = 100
        real(8)   ::   TOLERENCE_EPSON = 1.0E-14
  
    contains
    
        procedure  ::  CALCUATE_SEEDS_WEIGHTS   => CALCULATE_GLEGENDRE
        
    end type GAUSS_LEGENDRE_TYPE
    
    
    contains
    
    subroutine CALCULATE_GLEGENDRE( THIS  )
    
    class( GAUSS_LEGENDRE_TYPE) , intent(inout ) :: THIS
 
  
    !--------------LOCAL VARIABLES-----------------------!
    integer   ::   II, JJ, KK,ALLOCATE_ARRAY_FLAG
    integer   ::   MM
    real(8)   ::   Z_INPUT, Z1_INPUT, P_N, P_NMINUS1, P3,PP_DERIVATIVE
    
  
    
    MM = ( THIS%NUMBER_POINTS + 1 )/2  !---------NUMBER OF ROOTS TO COMPUTATE-------------!
    
 
    
    do II = 1, MM  !-----------------------LOOP OF HALP POINTS--------------!PI
        
        
        !-----------------INITIAL GUESS---------------------!
        Z_INPUT = DCOS( pi*(4.0D0*II - 1.0D0 )/(4.0D0*THIS%NUMBER_POINTS + 2.0D0 ) )
        
        
      do JJ = 1, THIS%NEWTON_ITEARTION_MAX
          
          !---------------COMPUTE THE LEGENDRE P_N(Z) USING RECURRENCE---------------!
          
          P_N = 1.0D0
          P_NMINUS1 = 0.0D0
          
          do KK = 1, THIS%NUMBER_POINTS
              
              P3 = P_NMINUS1
              P_NMINUS1 = P_N
              P_N = (  ( 2.0D0*KK - 1.0D0)*Z_INPUT*P_NMINUS1 - ( KK - 1.0D0)*P3 ) /KK
              
          end do
          ! Now p1 = P_n(z), p2 = P_{n-1}(z)
          ! Derivative P_n'(z)
          PP_DERIVATIVE = THIS%NUMBER_POINTS * (Z_INPUT*P_N - P_NMINUS1) / (Z_INPUT*Z_INPUT - 1.0d0)
          
          Z1_INPUT = Z_INPUT
          Z_INPUT = Z1_INPUT - P_N/PP_DERIVATIVE
          
          if( DABS(Z_INPUT - Z1_INPUT ) <= THIS%TOLERENCE_EPSON )then
              
              GOTO 115
          end if
      
      end do  !--------------NEWTON ALGORITHM-------------------!
      
115       THIS%SEEDS_POINT( II ) =  - Z_INPUT      
          THIS%SEEDS_POINT( THIS%NUMBER_POINTS + 1 - II ) =  + Z_INPUT  
          
          THIS%WEIGHTS_POINT( II ) = 2.0D0/( (1.0D0 - Z_INPUT*Z_INPUT )*PP_DERIVATIVE*PP_DERIVATIVE )
          THIS%WEIGHTS_POINT( THIS%NUMBER_POINTS + 1 - II ) = THIS%WEIGHTS_POINT( II )        
        
    end do  !-----------------------LOOP OF HALP POINTS--------------!
    
    
   THIS%SEEDS_FORMING_FLAG = .TRUE. 
    
  
    
    end subroutine CALCULATE_GLEGENDRE
    
    
    
    
    
    end module GAUSS_LEDENDRE_MODULE   
   
    
    
    
    
  