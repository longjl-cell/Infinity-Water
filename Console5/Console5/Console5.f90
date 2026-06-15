!=============================================================================!
!  fxy_integral.f90
!
!  Computes:
!    F(X,Y) = -2*exp(-Y) * INT_0^Y  e^t * (X^2 + t^2)^(-1/2)  dt
!             - pi*exp(-Y) * [Y0(X) + H0(X)]
!
!  where:
!    Y0(X) = Bessel function of the second kind, order zero
!    H0(X) = Struve function of order zero
!=============================================================================!
 

 
 

 
 
!=============================================================================!

 
 program  TEST_PROGRAM
 USE quadrature_INTEGRATION_mod
 use precision_mod
 use GAUSS_LEDENDRE_MODULE
 use special_functions_mod
 
 implicit none
 
 integer  ::  ALLCOAT_ARRAY_FLAG
 integer  ::  NUMBER_GAUSS_POINTS, II, JJ
 
 class(GAUSS_LEGENDRE_TYPE), allocatable :: GL
 class( SPECIAL_FUNCITON_TYPE), allocatable :: SPECIAL_FUNCTION
 
 real(DP), dimension(:), allocatable ::   X_VALUE
 real(DP), dimension(:), allocatable ::   Y_VALUE, OUTPUT_RESULT
 
 real(DP ) ::     SCALE1, H0, Y0
 
 
 
 
 
 
 
 allocate( GL, SPECIAL_FUNCTION,    stat = ALLCOAT_ARRAY_FLAG )
 allocate( X_VALUE(9), Y_VALUE( 8), OUTPUT_RESULT(8),  stat = ALLCOAT_ARRAY_FLAG )
 
 
 X_VALUE(:) = 0.0D0
 
 X_VALUE(1) = 0.10_DP
 X_VALUE(2) = 0.20_DP 
 X_VALUE(3) = 0.50_DP
 X_VALUE(4) = 1.00_DP
 X_VALUE(5) = 2.00_DP
 X_VALUE(6) = 5.00_DP
 X_VALUE(7) = 10.0_DP
 X_VALUE(8) = 20.0_DP
 X_VALUE(9) = 50.0_DP
 
 
 Y_VALUE(1) = 0.10_DP
 Y_VALUE(2) = 0.20_DP 
 Y_VALUE(3) = 0.50_DP
 Y_VALUE(4) = 1.00_DP
 Y_VALUE(5) = 2.00_DP
 Y_VALUE(6) = 5.00_DP
 Y_VALUE(7) = 10.0_DP
 Y_VALUE(8) = 20.0_DP
 
 
 
 
 NUMBER_GAUSS_POINTS = 30
 call GL%ALLCOATE_ARRAY( NUMBER_GAUSS_POINTS)
 
 CALL GL%CALCUATE_SEEDS_WEIGHTS( )

 
 !class( quadrature_INTEGRATION_mod) :: GL
 
 open(115, file = 'OUTPUT_RESULT.DAT' )
 
 do II = 1, 9
     
     do JJ = 1, 8
      
         
          SCALE1 = - PI*EXP( -Y_VALUE( JJ ) )
CALL SPECIAL_FUNCTION%BesselY0( X_VALUE( II ) , Y0)
call SPECIAL_FUNCTION%StruveH0( X_VALUE( II ) , H0)




       call integrate_kernel( X_VALUE( II ), Y_VALUE( JJ ), GL, OUTPUT_RESULT (JJ ) )   !------------------FIRST TERMS
       
       OUTPUT_RESULT (JJ )  = -2.0_DP*OUTPUT_RESULT ( JJ )  + SCALE1*( Y0 + H0)
       
       
     end do
     
     
     write( 115, 1150 ) ( OUTPUT_RESULT( JJ), JJ = 1, 8 )
     
 end do
 
 
1150 format( 1000(E30.15 )) 

  
 
  



    END program 
    
    
    
    
!=============================================================================!
!program fxy_main
!  use precision_mod
!  use special_functions_mod
!  use quadrature_mod
!  implicit none
! 
!  integer,  parameter :: N_GL = 64        ! GL points (high accuracy, no singularity)
! 
!  ! --- Grid matching the reference table ---
!  integer,  parameter :: NY = 8, NX = 9
!  real(dp), parameter :: Yvals(NY) = [0.1_dp, 0.2_dp, 0.5_dp, 1.0_dp, &
!                                       2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]
!  real(dp), parameter :: Xvals(NX) = [0.10_dp, 0.20_dp, 0.50_dp, 1.00_dp, &
!                                       2.00_dp, 5.00_dp, 10.00_dp, 20.00_dp, 50.00_dp]
! 
!  real(dp) :: X, Y, Ikernel, F_val
!  integer  :: ix, iy, iunit
! 
!  open(newunit=iunit, file='fxy_results.csv', status='replace', action='write')
!  write(iunit,'(A)', advance='no') 'Y'
!  do ix = 1, NX
!    write(iunit,'(A,F6.2)', advance='no') ',X=', Xvals(ix)
!  end do
!  write(iunit,*)
! 
!  write(*,'(A)') '======================================================================'
!  write(*,'(A)') '  F(X,Y) = -2*exp(-Y)*INT_0^Y e^t*(X^2+t^2)^(-1/2) dt'
!  write(*,'(A)') '           - pi*exp(-Y)*[Y0(X) + H0(X)]'
!  write(*,'(A)') '======================================================================'
! 
!  do iy = 1, NY
!    Y = Yvals(iy)
! 
!    write(iunit,'(F6.2)', advance='no') Y
!    write(*,'(A,F6.2)') '  Y = ', Y
! 
!    do ix = 1, NX
!      X = Xvals(ix)
! 
!      ! --- Term 1: -2*exp(-Y) * INT_0^Y e^t*(X^2+t^2)^(-1/2) dt ---
!      Ikernel = integrate_kernel(X, Y, N_GL)
! 
!      ! --- Term 2: -pi*exp(-Y)*[Y0(X) + H0(X)] ---
!      F_val = -2.0_dp * exp(-Y) * Ikernel &
!              - pi * exp(-Y) * (besy0(X) + struveh0(X))
! 
!      write(*,'(A,F6.2,A,ES22.15)') '    X = ', X, '   F = ', F_val
!      write(iunit,'(A,ES22.15)', advance='no') ',', F_val
!    end do
! 
!    write(*,*)
!    write(iunit,*)
!  end do
! 
!  close(iunit)
!  write(*,'(A)') 'Results written to fxy_results.csv'
! 
!end program fxy_main

