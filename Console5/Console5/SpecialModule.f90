!=============================================================================!
module special_functions_mod

  use precision_mod
  
  implicit none
 
  type SPECIAL_FUNCITON_TYPE
      
  integer  ::   NUMBER_TERM = 100000
integer  ::   NEWTON_ITERATION_NUMBER = 150
real(8)  ::   EPS_ERROR  =  1.0E-20
real(8)  ::   EPS_ERROR_NEWTON  =  1.0E-15
real(8)  ::   GAMA_CONSTANT = 0.57721566490153286060651209008240D0
real(8)  ::   STEP_SIZE = 0.00050D0   !----THIS PARAMETER IS IMPORTANT FOR THE ACCURACY----------!
  contains
  
       procedure   ::  StruveH0, StruveH1
       procedure   ::  BesselY0, BesselY1
       procedure   ::  BesselJ0, BesselJ1
       
       procedure  :: EXPONENT_INTEGRAL_E0_MODIFY

  end type  SPECIAL_FUNCITON_TYPE
  
  

    contains
  subroutine EXPONENT_INTEGRAL_E0_MODIFY( THIS, X_INPUT, X_OUTPUT )
 implicit none
class( SPECIAL_FUNCITON_TYPE ), intent(in) :: THIS
real(8), INTENT(in)     ::   X_INPUT
real(8), intent( out)   ::   X_OUTPUT

!-------------------LOCAL VARAIBLES-----------------------------!
integer   ::  KK
real(8)  :: X_INPUT_INVERSE, SUMM, TERM,ADD
real(8)  ::  KK_REAL, KK_REALP1, KK_REAL_INVERSE,KK_REALP1_INVERSE

    if (X_INPUT <= 0d0) then
       write(*,*)  'Ei(x) implemented for x > 0 only'
      PAUSE
    end if

 X_INPUT_INVERSE = 1.0D0/X_INPUT

    !------------------------------------------------
    ! Region 1: small/moderate x → series expansion
    !------------------------------------------------
    if (X_INPUT <= 40d0) then

       SUMM  = 0d0
       TERM = X_INPUT          ! term represents x^k / k!

       do KK = 1, THIS%NEWTON_ITERATION_NUMBER
           KK_REAL = real(KK)
           KK_REALP1 = KK_REAL + 1.0D0
           
           KK_REAL_INVERSE = 1.0D0/KK_REAL
           KK_REALP1_INVERSE = 1.0D0/KK_REALP1
           
          ADD = TERM*KK_REAL_INVERSE   ! x^k / (k*k!)
          SUMM = SUMM + ADD

        if (DABS( ADD ) < THIS%EPS_ERROR*DABS( SUMM + 1.0d0)) exit

          TERM = TERM * X_INPUT *KK_REALP1_INVERSE
       end do

       X_OUTPUT = THIS%GAMA_CONSTANT + DLOG( X_INPUT ) + SUMM

       X_OUTPUT = DEXP( - X_INPUT )*X_OUTPUT
    else
 !-----------------------------------------------------------
  ! Compute expEi_scaled(x) = exp(-x) * Ei(x) for x > 0
  ! using the large-x asymptotic expansion
  !
  ! exp(-x)Ei(x) ≈ (1/x) * sum_{m=0}^{M-1} m!/x^m
  !
  ! This is the safest way to handle very large x (avoids overflow).
  !-----------------------------------------------------------
  X_INPUT_INVERSE = 1.0d0/X_INPUT
  TERM = 1.0d0
  SUMM  = 1.0d0

  do KK = 0, 200
      
     TERM = TERM * real(KK+1 ) * X_INPUT_INVERSE      ! (m+1)! / x^(m+1) relative update
     SUMM  = SUMM + TERM
     if (abs(TERM) <= THIS%EPS_ERROR*abs(SUMM)) exit
     if (real(KK +2 ) > X_INPUT) exit            ! stop before divergence
  end do

  X_OUTPUT = X_INPUT_INVERSE * SUMM
 
    end if
    
   end subroutine EXPONENT_INTEGRAL_E0_MODIFY
    
 SUBROUTINE BesselY0( THIS, Y_INPUT, Y_OUTPUT )
implicit none
class( SPECIAL_FUNCITON_TYPE ), intent(in) :: THIS
real(DP), INTENT(in)  ::   Y_INPUT
real(DP), intent( out)   ::   Y_OUTPUT
  
! --- Local variables -------------------------------------
	REAL(DP)				::	f0,theta0
 	real(DP)             ::  TEMP_VALUE		

	if(Y_INPUT <= 3.0D0) then
        call THIS%BesselJ0(Y_INPUT, TEMP_VALUE )
		Y_OUTPUT	=	(2.0D0/pi)*dlog(Y_INPUT/2.0D0)*TEMP_VALUE	&	
						+0.367466907D0							&
						+0.605593797D0*(Y_INPUT/3.0D0)**2			&														
						-0.743505078D0*(Y_INPUT/3.0D0)**4			&														
						+0.253005481D0*(Y_INPUT/3.0D0)**6			&														
						-0.042619616D0*(Y_INPUT/3.0D0)**8			&														
						+0.004285691D0*(Y_INPUT/3.0D0)**10			&														
						-0.000250716D0*(Y_INPUT/3.0D0)**12														
	else
		f0			=	 0.79788454D0					&
						-0.00553897D0*(3.0D0/Y_INPUT)**2		&														
						+0.00099336D0*(3.0D0/Y_INPUT)**4		&														
						-0.00044346D0*(3.0D0/Y_INPUT)**6		&														
						+0.00020445D0*(3.0D0/Y_INPUT)**8		&														
						-0.00004959D0*(3.0D0/Y_INPUT)**10														

		theta0		=	 Y_INPUT		-	0.25D0*pi			&
						-0.04166592D0*(3.0D0/Y_INPUT)**1		&														
						+0.00239399D0*(3.0D0/Y_INPUT)**3		&														
						-0.00073984D0*(3.0D0/Y_INPUT)**5		&														
						+0.00031099D0*(3.0D0/Y_INPUT)**7		&														
						-0.00007605D0*(3.0D0/Y_INPUT)**9														

		Y_OUTPUT	=	f0*dsin(theta0)/dsqrt(Y_INPUT)																		

	end if
 
    end SUBROUTINE BesselY0
    
    
     SUBROUTINE BesselY1( THIS, Y_INPUT, Y_OUTPUT )
implicit none
class( SPECIAL_FUNCITON_TYPE ), intent(in) :: THIS
real(DP), INTENT(in)  ::   Y_INPUT
real(DP), intent( out)   ::   Y_OUTPUT

! --- Local variables -------------------------------------
	REAL(DP)				::	f1,theta1,TEMP_VALUE
 			

	if(Y_INPUT <= 3.0D0) then
         call THIS%BesselJ1(Y_INPUT, TEMP_VALUE )
		Y_OUTPUT	=	(2.0D0/pi)*(dlog(Y_INPUT/2.0D0)*TEMP_VALUE-1.0D0/Y_INPUT)	&	
						+0.07373571D0*(Y_INPUT/3.0D0)**1				&														
						+0.72276433D0*(Y_INPUT/3.0D0)**3				&														
						-0.43885620D0*(Y_INPUT/3.0D0)**5				&														
						+0.10418264D0*(Y_INPUT/3.0D0)**7				&														
						-0.01340825D0*(Y_INPUT/3.0D0)**9				&														
						+0.00094249D0*(Y_INPUT/3.0D0)**11														
	else
		f1			=	 0.79788459D0						&
						+0.01662008D0*(3.0D0/Y_INPUT)**2		&														
						-0.00187002D0*(3.0D0/Y_INPUT)**4		&														
						+0.00068519D0*(3.0D0/Y_INPUT)**6		&														
						-0.00029440D0*(3.0D0/Y_INPUT)**8		&														
						+0.00006952D0*(3.0D0/Y_INPUT)**10														

		theta1		=	 Y_INPUT		-	3.0D0*pi/4.0D0		&
						+0.12499895D0*(3.0D0/Y_INPUT)**1		&														
						-0.00605240D0*(3.0D0/Y_INPUT)**3		&														
						+0.00135825D0*(3.0D0/Y_INPUT)**5		&														
						-0.00049616D0*(3.0D0/Y_INPUT)**7		&														
						+0.00011531D0*(3.0D0/Y_INPUT)**9														

		Y_OUTPUT	=	f1*dsin(theta1)/dsqrt(Y_INPUT)																		

	end if
 
 end SUBROUTINE BesselY1
  
 
  SUBROUTINE StruveH1( THIS, Y_INPUT, Y_OUTPUT )
implicit none
class( SPECIAL_FUNCITON_TYPE ), intent(in) :: THIS
real(DP), INTENT(in)     ::   Y_INPUT
real(DP), intent( out)   ::   Y_OUTPUT
 
 

! --- Local variables -------------------------------------

	integer					::	ii
	REAL*8					::	P1,P2,P3
	REAL*8					::	P4,P5,P6
	REAL*8					::	a0,a1,a2,a3
	REAL*8					::	b1,b2,b3
	REAL*8					::	c1,c2,yy
 			

	if(Y_INPUT <= 3.0D0) then
        
		yy	=	(Y_INPUT/3.0D0)**2
        
		P1	=	+1.909859286D0
		P2	=	-1.145914713D0
		P3	=	+0.294656958D0
		P4	=	-0.042070508D0
 		P5	=	+0.003785727D0
		P6	=	-0.000207183D0
      
		Y_OUTPUT	=	(P1+(P2+(P3+(P4+(P5+P6*yy)*yy)*yy)*yy)*yy)*yy

    else
        
		yy	=	(3.0D0/Y_INPUT)**2

		a0	=	1.00000004D0
		a1	=	3.92205313D0
		a2	=	2.64893033D0
		a3	=	0.27450895D0

		b1	=	3.81095112D0
		b2	=	2.26216956D0
		b3	=	0.10885141D0

		c1	=	2.0D0*(a0	+	(a1+(a2+a3*yy)*yy)*yy) 
		c2	=	pi*(1.0D0	+	(b1+(b2+b3*yy)*yy)*yy) 

        call THIS%BesselY1(Y_INPUT, Y_OUTPUT )
		Y_OUTPUT	=	Y_OUTPUT + c1/c2	 

	end if

 
 end SUBROUTINE StruveH1
   
subroutine BesselJ0( THIS, Y_INPUT, Y_OUTPUT )
implicit none
class( SPECIAL_FUNCITON_TYPE ), intent(in) :: THIS
real(DP), INTENT(in)  ::   Y_INPUT
real(DP), intent( out)   ::   Y_OUTPUT
 
! --- Local variables -------------------------------------
	REAL( DP ) 	::	    yy,y2,f0,theta0			
	REAL( DP ) 	::	    P0,P1,P2,P3,P4,P5,P6
	REAL( DP ) 	::	    R0,R1,R2,R3,R4,R5
	REAL( DP ) 	::      S1,S2,S3,S4,S5
    
    
	data	P0,P1,P2,P3,P4,P5,P6		&
		/	+0.999999999D0,	-2.249999879D0,	+1.265623060D0,	&
			-0.316394552D0,	+0.044460948D0,	-0.003954479D0,	&
			+0.000212950D0	/
    
	data	R0,R1,R2,R3,R4,R5			&
		/	+0.79788454D0,	-0.00553897D0,	+0.00099336D0,	&
			-0.00044346D0,	+0.00020445D0,	-0.00004959D0	/
    
	data	S1,S2,S3,S4,S5				&
		/	-0.04166592D0,	+0.00239399D0,	-0.00073984D0,	&
			+0.00031099D0,	-0.00007605D0	/

Y_OUTPUT = 0.0D0

	if(Y_INPUT <= 3.0D0) then
        
		yy	=	(Y_INPUT/3.0D0)**2
		Y_OUTPUT	=	 P0+(P1+(P2+(P3+(P4+(P5+P6*yy)*yy)*yy)*yy)*yy)*yy
        
	else
		yy	=	3.0D0/Y_INPUT
		y2	=	yy**2
        
		f0			=	R0+(R1+(R2+(R3+(R4+R5*y2)*y2)*y2)*y2)*y2

		theta0		=	 Y_INPUT	-	0.25D0*pi		&
						+(S1+(S2+(S3+(S4+S5*y2)*y2)*y2)*y2)*yy

		Y_OUTPUT	=	f0*dcos(theta0)/dsqrt(Y_INPUT)
        
        
	end if

end subroutine BesselJ0 
 
subroutine BesselJ1( THIS, Y_INPUT, Y_OUTPUT )
implicit none
class( SPECIAL_FUNCITON_TYPE ), intent(in) :: THIS
real(DP), INTENT(in)  ::   Y_INPUT
real(DP), intent( out)   ::   Y_OUTPUT
 
! --- Local variables -------------------------------------
	REAL( 8 )					::	yy,y2,f1,theta1
	REAL( 8 )					::	P0,P1,P2,P3,P4,P5,P6
	REAL( 8 )					::	R0,R1,R2,R3,R4,R5
	REAL( 8 )					::  S1,S2,S3,S4,S5
    
	data	P0,P1,P2,P3,P4,P5,P6		&
		/	+0.500000000D0,	-0.562499992D0,	+0.210937377D0,	&
			-0.039550040D0,	+0.004447331D0,	-0.000330547D0,	&
			+0.000015525D0	/
    
	data	R0,R1,R2,R3,R4,R5			&
		/	+0.79788459D0,	+0.01662008D0,	-0.00187002D0,	&
			+0.00068519D0,	-0.00029440D0,	+0.00006952D0	/
    
	data	S1,S2,S3,S4,S5				&
		/	+0.12499895D0,	-0.00605240D0,	+0.00135825D0,	&
			-0.00049616D0,	+0.00011531D0	/

	if(Y_INPUT <= 3.0D0) then
		yy	=	Y_INPUT/3.0D0
		y2	=	yy*yy
		Y_OUTPUT	=	P0+(P1+(P2+(P3+(P4+(P5+P6*y2)*y2)*y2)*y2)*y2)*y2

		Y_OUTPUT	=	Y_OUTPUT*Y_INPUT
	else
		yy	=	3.0D0/Y_INPUT
		y2	=	yy*yy
		f1			=	R0+(R1+(R2+(R3+(R4+R5*y2)*y2)*y2)*y2)*y2

		theta1		=	 Y_INPUT	-	0.75D0*pi		&
						+(S1+(S2+(S3+(S4+S5*y2)*y2)*y2)*y2)*yy

		Y_OUTPUT	=	f1*dcos(theta1)/dsqrt(Y_INPUT)
	end if



end subroutine BesselJ1  


subroutine StruveH0( THIS, Y_INPUT, Y_OUTPUT )
implicit none
class( SPECIAL_FUNCITON_TYPE ), intent(in) :: THIS
real(DP), INTENT(in)  ::   Y_INPUT
real(DP), intent( out)   ::   Y_OUTPUT
 
! --- Local variables -------------------------------------

	integer					::	ii
	REAL( DP)					::	P0,P1,P2
	REAL( DP)					::	P3,P4,P5
	REAL( DP)					::	a0,a1,a2,a3
	REAL( DP)					::	b1,b2,b3
	REAL( DP)					::	c1,c2
	REAL( DP)					::	yy 

	if(Y_INPUT <= 3.0D0) then
        
		yy	=	(Y_INPUT/3.0D0)**2
        
		P0	=	+1.909859164D0
		P1	=	-1.909855001D0
		P2	=	+0.687514637D0
		P3	=	-0.126164557D0
 		P4	=	+0.013828813D0
		P5	=	-0.000876918D0
        
		Y_OUTPUT	=	P0+(P1+(P2+(P3+(P4+P5*yy)*yy)*yy)*yy)*yy
        
		Y_OUTPUT	=	Y_OUTPUT*(Y_INPUT/3.0D0)
        
    else
        
		yy	=	(3.0D0/Y_INPUT)**2

		a0	=	0.99999906D0
		a1	=	4.77228920D0
		a2	=	3.85542044D0
		a3	=	0.32303607D0

		b1	=	4.88331068D0
		b2	=	4.28957333D0
		b3	=	0.52120508D0

		c1	=	2.0D0*(a0	+	(a1+(a2+a3*yy)*yy)*yy) 
		c2	=	pi*Y_INPUT*(1.0D0+	(b1+(b2+b3*yy)*yy)*yy) 
							
        CALL THIS%BesselY0(Y_INPUT, Y_OUTPUT)		
        
		Y_OUTPUT	=	Y_OUTPUT + c1/c2													
																				
	end if
 
end subroutine  StruveH0
 
end module special_functions_mod