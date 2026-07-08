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
 use GREEN_TRANSFORM_MOD
 
 implicit none
 
 integer  ::  ALLCOAT_ARRAY_FLAG
 integer  ::  NUMBER_GAUSS_POINTS, II, JJ
 integer  ::  METHODOLOGY_FLAG
 
 class(GAUSS_LEGENDRE_TYPE), allocatable :: GL
  class(GAUSS_LAGUERRE_TYPE), allocatable :: GLA
 class( SPECIAL_FUNCITON_TYPE), allocatable :: SPECIAL_FUNCTION
 
 real(DP), dimension(:), allocatable ::   X_VALUE, Y_VALUE, OUTPUT_RESULT,GREEN_FUNCTION
 
 real(DP ) ::     SCALE1, H0, Y0
 
 
 real(dp) :: KERNEL_OUT(6)
 
 real(DP) :: K0_DEMO, X_PHYS, Y_PHYS, Z_PHYS, ZETA_DEMO, ETA_DEMO, XI_DEMO
 real(DP) :: G_DEMO, GX_DEMO, GY_DEMO, GZ_DEMO, GZETA_DEMO, GETA_DEMO, GXI_DEMO
 real(DP) :: GXX_DEMO, GYY_DEMO, GXY_DEMO, GZZ_DEMO, GXZ_DEMO, GYZ_DEMO
 
 
 
 
 allocate( GL, SPECIAL_FUNCTION, GLA,   stat = ALLCOAT_ARRAY_FLAG )
 allocate( X_VALUE(9), Y_VALUE( 8), OUTPUT_RESULT(8), GREEN_FUNCTION(6), stat = ALLCOAT_ARRAY_FLAG )
 
 
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
 
 
 
 
 NUMBER_GAUSS_POINTS = 100
 call GL%ALLCOATE_ARRAY( NUMBER_GAUSS_POINTS)
 
 CALL GL%CALCUATE_SEEDS_WEIGHTS( )
 
 
CALL GLA%ALLCOATE_ARRAY( NUMBER_GAUSS_POINTS)

  CALL GLA%CALCUATE_SEEDS_WEIGHTS( )
 !class( quadrature_INTEGRATION_mod) :: GL
 
 open(115, file = 'OUTPUT_RESULT.DAT' )
 
 METHODOLOGY_FLAG = 3
 
 
 select case( METHODOLOGY_FLAG)
 case(1)  !---------------------------------------------method 1----------------------------!
     open(117, file = 'METHOD1_RESULTS.csv')
    
     write(115, '(A)') '=== METHOD 1: Gauss-Legendre ==='
      write(115, '(A30, 6A30)') 'X', 'Y', 'GREEN', 'GREEN_DX', &
                                'GREEN_DY', 'GREEN_DXDX', &
                                'GREEN_DXDY', 'GREEN_DYDY'
      
       write(117, '(A)') 'X,Y,GREEN,GREEN_DX,GREEN_DY,GREEN_DXDX,GREEN_DXDY,GREEN_DYDY'


      do II = 1, 9
          do JJ = 1, 8

              call integrate_kernel( X_VALUE(II), Y_VALUE(JJ), &
                                      GL, SPECIAL_FUNCTION, KERNEL_OUT )

              write(115, 1150) X_VALUE(II), Y_VALUE(JJ), &
                               KERNEL_OUT(1), KERNEL_OUT(2), &
                               KERNEL_OUT(3), KERNEL_OUT(4), &
                               KERNEL_OUT(5), KERNEL_OUT(6)
              
            write(117, '(F10.4, A, F10.4, A, 6(E30.15, A))') &
                  X_VALUE(II), ',', Y_VALUE(JJ), ',', &
                  KERNEL_OUT(1), ',', KERNEL_OUT(2), ',', &
                  KERNEL_OUT(3), ',', KERNEL_OUT(4), ',', &
                  KERNEL_OUT(5), ',', KERNEL_OUT(6), ','
              
        
          end do
      end do
      
      close(117)
 
case(2)
    ! Write CSV header
    open(116, file = 'GREEN_RESULTS.csv')
    write(115, '(A)') '=== METHOD 2: Gauss-Laguerre Green function ==='
    write(116, '(A)') 'X,Y,GREEN,GREEN_DX,GREEN_DY,GREEN_DXDX,GREEN_DXDY,GREEN_DYDY'

    ! Write DAT header
    write(115, '(A30, 6A30)') 'X', 'Y', 'GREEN', 'GREEN_DX', &
                               'GREEN_DY', 'GREEN_DXDX', &
                               'GREEN_DXDY', 'GREEN_DYDY'

    do II = 1, 9
        do JJ = 1, 8

            call integral_kernel_laguerre( X_VALUE(II), Y_VALUE(JJ), &
                                           GLA, SPECIAL_FUNCTION, GREEN_FUNCTION)

            ! Write one row per (X,Y) pair to DAT file
            write(115, 1150) X_VALUE(II), Y_VALUE(JJ), &
                             GREEN_FUNCTION(1), GREEN_FUNCTION(2), &
                             GREEN_FUNCTION(3), GREEN_FUNCTION(4), &
                             GREEN_FUNCTION(5), GREEN_FUNCTION(6)

            ! Write same row to CSV
            write(116, '(F10.4, A, F10.4, A, 6(E30.15, A))') &
                X_VALUE(II), ',', Y_VALUE(JJ), ',', &
                GREEN_FUNCTION(1), ',', GREEN_FUNCTION(2), ',', &
                GREEN_FUNCTION(3), ',', GREEN_FUNCTION(4), ',', &
                GREEN_FUNCTION(5), ',', GREEN_FUNCTION(6), ','

        end do
    end do

    close(116)
     
     
    
         case(3)  !------------------------------ physical-space transform demo -----------------!
     ! This case does NOT sweep new physics - it reuses the same F(X,Y) grid as
     ! case(1), but shows where TRANSFORM_GREEN_DERIVATIVES belongs: it needs a
     ! REAL geometry (field point x,y,z and source/image point zeta,eta,xi), not
     ! just X,Y in isolation, because X,Y alone don't uniquely determine those.
     !
     ! Demo geometry: source/image at the origin (zeta=eta=xi=0), field point
     ! restricted to the x-axis (y=0) so that X = k0*x and Y = -k0*z invert
     ! cleanly for verification.
 
     open(118, file = 'METHOD3_TRANSFORM_RESULTS.csv')
     write(118, '(A)') 'X,Y,K0,x,z,G,GX,GY,GZ,GXX,GYY,GXY,GZZ,GXZ,GYZ'
 
     K0_DEMO   = 1.0_DP     ! pick your actual wavenumber here
     ZETA_DEMO = 0.0_DP
     ETA_DEMO  = 0.0_DP
     XI_DEMO   = 0.0_DP
 
     do II = 1, 9
         do JJ = 1, 8
 
             ! 1) evaluate F(X,Y) and its X,Y-derivatives exactly as in case(1)
             call integrate_kernel( X_VALUE(II), Y_VALUE(JJ), &
                                     GL, SPECIAL_FUNCTION, KERNEL_OUT )
 
             ! 2) recover the physical field point consistent with this (X,Y)
             !    under the demo geometry (y=0, source at origin)
             X_PHYS = X_VALUE(II) / K0_DEMO
             Y_PHYS = 0.0_DP
             Z_PHYS = -Y_VALUE(JJ) / K0_DEMO
 
             ! 3) apply the chain-rule transform to physical derivatives
             call TRANSFORM_GREEN_DERIVATIVES( X_PHYS, Y_PHYS, Z_PHYS, &
                                                ZETA_DEMO, ETA_DEMO, XI_DEMO, K0_DEMO, &
                                                KERNEL_OUT(1), KERNEL_OUT(2), KERNEL_OUT(3), &
                                                KERNEL_OUT(4), KERNEL_OUT(5), KERNEL_OUT(6), &
                                                G_DEMO, GX_DEMO, GY_DEMO, GZ_DEMO, &
                                                GZETA_DEMO, GETA_DEMO, GXI_DEMO, &
                                                GXX_DEMO, GYY_DEMO, GXY_DEMO, &
                                                GZZ_DEMO, GXZ_DEMO, GYZ_DEMO )
 
             write(118, '(2(F10.4,A), F6.2, A, 2(F10.4,A), 9(E30.15,A))') &
                 X_VALUE(II), ',', Y_VALUE(JJ), ',', K0_DEMO, ',', &
                 X_PHYS, ',', Z_PHYS, ',', &
                 G_DEMO, ',', GX_DEMO, ',', GY_DEMO, ',', GZ_DEMO, ',', &
                 GXX_DEMO, ',', GYY_DEMO, ',', GXY_DEMO, ',', &
                 GZZ_DEMO, ',', GXZ_DEMO, ',', GYZ_DEMO, ','
 
         end do
     end do
 
     close(118)
     
 
     
     
     
 case default
     write(*,*) 'there is no this option, please contact the author'
     stop
     
     end select 
 
 close(115)
 

 
 
1150 format( 8(E30.15 )) 

  
 
  



    END program 
    
    
    
    
!!=============================================================================!
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

