
    
  module quadrature_INTEGRATION_mod
  use precision_mod
  use GAUSS_LEDENDRE_MODULE
  implicit none
 
  public ::  integrate_kernel
  
  
  
contains
 
 FUNCTION FUNCTION_INTEGRAND(   T_BAR , X_INPUT, Y_INPUT ) RESULT( FUNCTION_VALUE )
    real(DP), intent( in  )     ::   T_BAR ,X_INPUT,Y_INPUT
    REAL(8)   ::  FUNCTION_VALUE
    !----------------------LOCAL VARIABLES-----------------!

    
 FUNCTION_VALUE = DEXP(T_BAR - Y_INPUT )/DSQRT(  X_INPUT*X_INPUT +T_BAR*T_BAR  )
 
 
 
 end function FUNCTION_INTEGRAND
 
  !---------------------------------------------------------------------------!
  ! Integrate  e^t * (X^2 + t^2)^(-1/2)  over [0, Y]
  ! using Gauss-Legendre quadrature
  !---------------------------------------------------------------------------!
  subroutine  integrate_kernel(X, Y, GAUSS_LEGENDER, OUTPUT_RESULT)  
    real(dp), intent(in) :: X, Y
    class(GAUSS_LEGENDRE_TYPE),  intent(in) :: GAUSS_LEGENDER
    real(DP ), intent( out )  :: OUTPUT_RESULT
    
    
    real(dp)             :: Ival
 
 
    real(DP) ::  A_VALUE, B_VALUE, SCALE_FACTOR,T_INPUT,FUNCTION_VALUE
    
    
    integer  :: I
 
  
    
     A_VALUE = 0.0_dp
     B_VALUE = Y
     
     
    Ival = 0.0_dp
    
    do I = 1, GAUSS_LEGENDER%NUMBER_POINTS   !----------------------LOOP OF GAUSS INTEGRATION-------------------------------!
        
         T_INPUT = 0.50D0*( B_VALUE + A_VALUE ) + 0.50D0*(  B_VALUE - A_VALUE )*GAUSS_LEGENDER%SEEDS_POINT( I )     
         
         
        SCALE_FACTOR =  0.50D0*(  B_VALUE  - A_VALUE )
        
        
        
       FUNCTION_VALUE = FUNCTION_INTEGRAND(  T_INPUT, X, Y )
      
      Ival = Ival + SCALE_FACTOR*GAUSS_LEGENDER%WEIGHTS_POINT( I )*FUNCTION_VALUE
      
      
    end do   !----------------------LOOP OF GAUSS INTEGRATION-------------------------------!
    
    
   
 OUTPUT_RESULT = Ival
 
    
  end subroutine integrate_kernel
  
  
 
   
 
 
end module quadrature_INTEGRATION_mod