
    
  module quadrature_INTEGRATION_mod
  use precision_mod
  use GAUSS_LEDENDRE_MODULE
  use GAUSS_LAGUREE_MODULE
  use special_functions_mod
  implicit none
 
  public ::  integrate_kernel 
  public ::  integral_kernel_laguerre
  
  
  
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
        
         T_INPUT = 0.50D0*( B_VALUE + A_VALUE ) + 0.50D0*(  B_VALUE - A_VALUE )*GAUSS_LEGENDER%SEEDS_POINT( I )  !-------------CONVERT VARIABINE FROM (0,Y) INTO (-1, 1)-----------------!   
         
         
        SCALE_FACTOR =  0.50D0*(  B_VALUE  - A_VALUE )
        
        
        
       FUNCTION_VALUE = FUNCTION_INTEGRAND(  T_INPUT, X, Y )
      
      Ival = Ival + SCALE_FACTOR*GAUSS_LEGENDER%WEIGHTS_POINT( I )*FUNCTION_VALUE
      
      
    end do   !----------------------LOOP OF GAUSS INTEGRATION-------------------------------!
    
    
   
 OUTPUT_RESULT = Ival
 
    
  end subroutine integrate_kernel
  
  
  subroutine integral_kernel_laguerre(X, Y,   GAUSS_LAGUERRE,SPECIAL_FUNCTION, OUTPUT_RESULT)  
  real(dp), intent(in) :: X, Y 
  class(GAUSS_LAGUERRE_TYPE),  intent(in) :: GAUSS_LAGUERRE
  class(SPECIAL_FUNCITON_TYPE), intent(in) ::  SPECIAL_FUNCTION
  real(dp ), dimension(6), intent( out )  :: OUTPUT_RESULT 
  
  integer :: II
  REAL(8) :: K_INPUT,K0_BAR, X_BAR, Y_BAR
  real(8) ::GREEN, GREEN_DX, GREEN_DY, GREEN_DXDX, GREEN_DXDY, GREEN_DYDY, SCALE_FACTOR, FUNCTION_VALUE,PARA_ALPHA,PARA_ALPHA_INVERSE, K0_BAR_INPUT, TEMP_VALUE1,EXPONENT_INTEGRAL_E0
  
  !----------------------NONDIMENSIONLIZED DATA---------------------------------------!
  X_BAR = X
  Y_BAR = Y

  
  K0_BAR = 1.0D0
  
  PARA_ALPHA =   Y_BAR
  PARA_ALPHA_INVERSE = 1.0/PARA_ALPHA
  GREEN = 0.0D0
  GREEN_DX = 0.0D0
  GREEN_DY = 0.0D0
  GREEN_DXDX = 0.0D0
  GREEN_DXDY = 0.0D0
  GREEN_DYDY = 0.0D0
  
  SCALE_FACTOR = 2.0*PARA_ALPHA_INVERSE
  
  
  !------------------------------------CALCUALTE I1------------------------------!
  do II = 1, GAUSS_LAGUERRE%NUMBER_POINTS !---------------------NUMBER OF GUASS INTEGRAITON POINTS-----------------------------------!
      
        K_INPUT =  PARA_ALPHA_INVERSE*GAUSS_LAGUERRE%SEEDS_POINT( II )   !-------------K/ALPHA--------------!
       
        
       FUNCTION_VALUE =  FUNCTION_COMBINE( K_INPUT, K0_BAR,X_BAR, SPECIAL_FUNCTION, 1 )   !----------------GREEN FUNCTION------------!
       GREEN = GREEN + SCALE_FACTOR*GAUSS_LAGUERRE%WEIGHTS_POINT( II )*FUNCTION_VALUE
       
     
       FUNCTION_VALUE =  FUNCTION_COMBINE( K_INPUT, K0_BAR,X_BAR, SPECIAL_FUNCTION, 2 )   !----------------GREEN FUNCTION------------!
       GREEN_DX = GREEN_DX + SCALE_FACTOR*GAUSS_LAGUERRE%WEIGHTS_POINT( II )*FUNCTION_VALUE
       
       FUNCTION_VALUE =  FUNCTION_COMBINE( K_INPUT, K0_BAR,X_BAR, SPECIAL_FUNCTION, 3 )   !----------------GREEN FUNCTION------------!
       GREEN_DY = GREEN_DY + SCALE_FACTOR*GAUSS_LAGUERRE%WEIGHTS_POINT( II )*FUNCTION_VALUE  
       
       
      FUNCTION_VALUE =  FUNCTION_COMBINE( K_INPUT, K0_BAR,X_BAR, SPECIAL_FUNCTION, 4 )   !----------------GREEN FUNCTION------------!
       GREEN_DXDX = GREEN_DXDX + SCALE_FACTOR*GAUSS_LAGUERRE%WEIGHTS_POINT( II )*FUNCTION_VALUE  
       
           FUNCTION_VALUE =  FUNCTION_COMBINE( K_INPUT, K0_BAR,X_BAR, SPECIAL_FUNCTION, 5 )   !----------------GREEN FUNCTION------------!
       GREEN_DXDY = GREEN_DXDY + SCALE_FACTOR*GAUSS_LAGUERRE%WEIGHTS_POINT( II )*FUNCTION_VALUE  
       
             FUNCTION_VALUE =  FUNCTION_COMBINE( K_INPUT, K0_BAR,X_BAR, SPECIAL_FUNCTION, 6 )   !----------------GREEN FUNCTION------------!
       GREEN_DYDY = GREEN_DYDY + SCALE_FACTOR*GAUSS_LAGUERRE%WEIGHTS_POINT( II )*FUNCTION_VALUE  
      
  end do  !---------------------NUMBER OF GUASS INTEGRAITON POINTS-----------------------------------!
  
  
  
  
     K0_BAR = K0_BAR   !-----------------------k0--------------------!  
   K0_BAR_INPUT = PARA_ALPHA*K0_BAR         !----------------------alpha------------------!
   
  call  SPECIAL_FUNCTION%EXPONENT_INTEGRAL_E0_MODIFY( K0_BAR_INPUT, EXPONENT_INTEGRAL_E0)  !-----------EXPONENTENTAL INTEGRAL----(INPUT AS k0)----------J1---------!
 
   TEMP_VALUE1 = FUNCTION_COMBINATION_NEWFORMULATIONPART2(  K0_BAR ,X_BAR,SPECIAL_FUNCTION, 1, PARA_ALPHA )                 !----------------GREEN-----------------------!
  
   GREEN = GREEN - 2.0D0*TEMP_VALUE1*EXPONENT_INTEGRAL_E0
   
   TEMP_VALUE1 = FUNCTION_COMBINATION_NEWFORMULATIONPART2(  K0_BAR ,X_BAR,SPECIAL_FUNCTION, 2, PARA_ALPHA )                 !----------------GREEN-----------------------!
   
    GREEN_DX = GREEN_DX - 2.0D0*TEMP_VALUE1*EXPONENT_INTEGRAL_E0
   

    TEMP_VALUE1 = FUNCTION_COMBINATION_NEWFORMULATIONPART2(  K0_BAR ,X_BAR, SPECIAL_FUNCTION,3, PARA_ALPHA )                 !----------------GREEN-----------------------!
   
    GREEN_DY = GREEN_DY- 2.0D0*TEMP_VALUE1*EXPONENT_INTEGRAL_E0
    
      TEMP_VALUE1 = FUNCTION_COMBINATION_NEWFORMULATIONPART2(  K0_BAR ,X_BAR, SPECIAL_FUNCTION,4, PARA_ALPHA )                 !----------------GREEN-----------------------!
   
    GREEN_DXDX = GREEN_DXDX- 2.0D0*TEMP_VALUE1*EXPONENT_INTEGRAL_E0  
    
    
      TEMP_VALUE1 = FUNCTION_COMBINATION_NEWFORMULATIONPART2(  K0_BAR ,X_BAR, SPECIAL_FUNCTION,5, PARA_ALPHA )                 !----------------GREEN-----------------------!
   
    GREEN_DXDY = GREEN_DXDY- 2.0D0*TEMP_VALUE1*EXPONENT_INTEGRAL_E0      
    
    
    
  TEMP_VALUE1 = FUNCTION_COMBINATION_NEWFORMULATIONPART2(  K0_BAR ,X_BAR, SPECIAL_FUNCTION,6, PARA_ALPHA )                 !----------------GREEN-----------------------!
   
    GREEN_DYDY = GREEN_DYDY- 2.0D0*TEMP_VALUE1*EXPONENT_INTEGRAL_E0  
   
 
  
  OUTPUT_RESULT(1) = GREEN
  OUTPUT_RESULT(2) = GREEN_DX
  OUTPUT_RESULT(3) = GREEN_DY
  OUTPUT_RESULT(4) = GREEN_DXDX
  OUTPUT_RESULT(5) = GREEN_DXDY
  OUTPUT_RESULT(6) = GREEN_DYDY
       
  
   
  
  
  
  end subroutine integral_kernel_laguerre
  
 
   
 function  FUNCTION_COMBINATION_NEWFORMULATIONPART2(  K0_BAR, R_INPUT, SPECIAL_FUNCTION, SWITH_FLAG , ALPHA_PARA )  result(FUNCTION_VALUE)
 real(8), intent( in ) ::  K0_BAR, R_INPUT
  class(SPECIAL_FUNCITON_TYPE), intent(in) ::  SPECIAL_FUNCTION
 integer, intent( in ) :: SWITH_FLAG
 real(8), intent( in )  :: ALPHA_PARA
 
 
 real(8)  ::  FUNCTION_VALUE
 
 
 !------------------------LOCAL VARIABLES------------------------!
 
 real(8) ::    K0R_INPUT,  J0_VALUE_K0, J1_VALUE_K0
 real(8) :: F_FUNCTION_1,  DENOMINATOR_G1 
 
 
 FUNCTION_VALUE = 0.0D0
 
 K0R_INPUT  = K0_BAR*R_INPUT
 
   call SPECIAL_FUNCTION%BesselJ0( K0R_INPUT, J0_VALUE_K0 )
   call SPECIAL_FUNCTION%BesselJ0( K0R_INPUT, J1_VALUE_K0 ) 
   
   F_FUNCTION_1 = 0.0

 select case( SWITH_FLAG)
 case(1)  !-------------F_FUCNTION-----------------------!
      F_FUNCTION_1 =  F_FUNCTION(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )
 
 case(2)  !------------DF/DX-----------------------!
     
      F_FUNCTION_1 =  F_FUNCTION_DX(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )
  
  case(3)  !------------DF/DY-----------------------!
     
      F_FUNCTION_1 =  F_FUNCTION_DY(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )
      
      
 case(4)  !------------D^2F/DX^2-----------------------!
      F_FUNCTION_1 =  F_FUNCTION_DXDX(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )

      
 case(5)  !------------D^2F/DXDY-----------------------!
     
      F_FUNCTION_1 =  F_FUNCTION_DXDY(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )
   case(6)  !------------D^2F/DX^2-----------------------!
       
      F_FUNCTION_1 =  F_FUNCTION_DYDY(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )
      
 end select
 
 

 DENOMINATOR_G1 = G_FUNCTION_1STDERIVATIVE(K0_BAR)
 
 
 FUNCTION_VALUE = F_FUNCTION_1/DENOMINATOR_G1
 
 
 
 end function FUNCTION_COMBINATION_NEWFORMULATIONPART2
 
 
 
 
 
 
 
 
 
 function  FUNCTION_COMBINE( K_BAR, K0_BAR, R_INPUT, SPECIAL_FUNCTION, SWITH_FLAG )  result(FUNCTION_VALUE)
 real(8), intent( in ) :: K_BAR, K0_BAR, R_INPUT
  class(SPECIAL_FUNCITON_TYPE), intent(in) ::  SPECIAL_FUNCTION
 integer, intent( in ) :: SWITH_FLAG
 real(8)  ::  FUNCTION_VALUE
 
 
 !------------------------LOCAL VARIABLES------------------------!
 
 real(8) ::  KR_INPUT, K0R_INPUT, J0_VALUE, J1_VALUE, J0_VALUE_K0, J1_VALUE_K0
 real(8) :: F_FUNCTION_1, F_FUNCTION_2, TEMP_VALUE, DENOMINATOR_G1, DENOMINATOR_G2
 
 
 FUNCTION_VALUE = 0.0D0
 
 
 KR_INPUT = K_BAR*R_INPUT
 
   call SPECIAL_FUNCTION%BesselJ0( KR_INPUT, J0_VALUE )
   call SPECIAL_FUNCTION%BesselJ0( KR_INPUT, J1_VALUE )
 
 K0R_INPUT  = K0_BAR*R_INPUT
 
   call SPECIAL_FUNCTION%BesselJ0( K0R_INPUT, J0_VALUE_K0 )
   call SPECIAL_FUNCTION%BesselJ0( K0R_INPUT, J1_VALUE_K0 ) 
   
   F_FUNCTION_1 = 0.0
   F_FUNCTION_2 = 0.0
   TEMP_VALUE = K_BAR - K0_BAR
   
 select case( SWITH_FLAG)
 case(1)  !-------------F_FUCNTION-----------------------!
      F_FUNCTION_1 =  F_FUNCTION(K_BAR, R_INPUT, J0_VALUE, J1_VALUE )
      F_FUNCTION_2 =  F_FUNCTION(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )    
 case(2)  !------------DF/DX-----------------------!
     
      F_FUNCTION_1 =  F_FUNCTION_DX(K_BAR, R_INPUT, J0_VALUE, J1_VALUE )
      F_FUNCTION_2 =  F_FUNCTION_DX(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )    
      
 case(4)  !------------D^2F/DX^2-----------------------!
      F_FUNCTION_1 =  F_FUNCTION_DXDX(K_BAR, R_INPUT, J0_VALUE, J1_VALUE )
      F_FUNCTION_2 =  F_FUNCTION_DXDX(K0_BAR, R_INPUT, J0_VALUE_K0, J1_VALUE_K0 )    
 end select
 
 
 
 DENOMINATOR_G1 = G_FUNCTION(K_BAR )
 DENOMINATOR_G2 = G_FUNCTION_1STDERIVATIVE(K0_BAR)
 
 TEMP_VALUE = 1.0/TEMP_VALUE
 F_FUNCTION_1 = F_FUNCTION_1/DENOMINATOR_G1
 F_FUNCTION_2 = TEMP_VALUE*F_FUNCTION_2/DENOMINATOR_G2
 
 FUNCTION_VALUE = F_FUNCTION_1 - F_FUNCTION_2
 
 
 
 end function FUNCTION_COMBINE
 
  function  F_FUNCTION(K_BAR, Y_INPUT, J0_VALUE, J1_VALUE )  RESULT(FUNCTION_VALUE)
 real(8) , INTENT( in) :: K_BAR, Y_INPUT, J0_VALUE, J1_VALUE
 real(8) :: FUNCTION_VALUE
 

 
 FUNCTION_VALUE = J0_VALUE  
 
 
 end function F_FUNCTION
 
  function  F_FUNCTION_DX( K_BAR, Y_INPUT, J0_VALUE, J1_VALUE )  RESULT(FUNCTION_VALUE)
 real(8),INTENT( in) :: K_BAR, Y_INPUT, J0_VALUE, J1_VALUE
 real(8) :: FUNCTION_VALUE
 

 
 FUNCTION_VALUE = - K_BAR*J1_VALUE  
 
 
  end function F_FUNCTION_DX
  
  
  function  F_FUNCTION_DY( K_BAR, Y_INPUT, J0_VALUE, J1_VALUE )  RESULT(FUNCTION_VALUE)
 real(8),INTENT( in) :: K_BAR, Y_INPUT, J0_VALUE, J1_VALUE
 real(8) :: FUNCTION_VALUE
 

 
 FUNCTION_VALUE = - K_BAR*J0_VALUE  
 
 
 end function F_FUNCTION_DY 
 
   function  F_FUNCTION_DXDX(K_BAR, Y_INPUT, J0_VALUE, J1_VALUE )  RESULT(FUNCTION_VALUE)
 real(8) ,INTENT( in) :: K_BAR, Y_INPUT, J0_VALUE, J1_VALUE
 real(8)  :: FUNCTION_VALUE
 
 !------------REAL
 real(8)::  KR_BAR_INVERSE

 KR_BAR_INVERSE = 1.0/(K_BAR*Y_INPUT)
 FUNCTION_VALUE = - K_BAR*K_BAR*(J0_VALUE  - KR_BAR_INVERSE*J1_VALUE )
 
 
   end function F_FUNCTION_DXDX
   
    function  F_FUNCTION_DXDY(K_BAR, Y_INPUT, J0_VALUE, J1_VALUE )  RESULT(FUNCTION_VALUE)
 real(8) ,INTENT( in) :: K_BAR, Y_INPUT, J0_VALUE, J1_VALUE
 real(8)  :: FUNCTION_VALUE
 
 !------------REAL
 real(8)::  KR_BAR_INVERSE


 FUNCTION_VALUE =  K_BAR*K_BAR*J1_VALUE 
 
 
    end function F_FUNCTION_DXDY
   
    
        function  F_FUNCTION_DYDY(K_BAR, Y_INPUT, J0_VALUE, J1_VALUE )  RESULT(FUNCTION_VALUE)
 real(8) ,INTENT( in) :: K_BAR, Y_INPUT, J0_VALUE, J1_VALUE
 real(8)  :: FUNCTION_VALUE
 
 !------------REAL
 real(8)::  KR_BAR_INVERSE


 FUNCTION_VALUE =  K_BAR*K_BAR*J0_VALUE 
 
 
   end function F_FUNCTION_DYDY
   
     function  G_FUNCTION(K_BAR ) result(FUNCTION_VALUE)
 real(8), INTENT( in) :: K_BAR
 real(8)  :: FUNCTION_VALUE
 

 
 FUNCTION_VALUE = K_BAR - 1.0  
 
 
     end function G_FUNCTION
     
          function  G_FUNCTION_1STDERIVATIVE(K_BAR )  RESULT(FUNCTION_VALUE)
 real(8),INTENT( in) :: K_BAR
 real(8) :: FUNCTION_VALUE
 

 
 FUNCTION_VALUE = 1.0
 
 
 end function G_FUNCTION_1STDERIVATIVE
 
end module quadrature_INTEGRATION_mod