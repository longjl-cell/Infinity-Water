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
 
module precision_mod
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793_dp
end module precision_mod
 
 
!=============================================================================!
module special_functions_mod
  use precision_mod
  implicit none
 
contains
 
  !---------------------------------------------------------------------------!
  ! Y0(x) : Bessel function of the second kind, order zero
  ! Abramowitz & Stegun 9.4.2 / 9.4.4
  !---------------------------------------------------------------------------!
  function besy0(x) result(res)
    real(dp), intent(in) :: x
    real(dp)             :: res
    real(dp)             :: z, xx
 
    real(dp), parameter :: rp1(7) = [ &
        -2957821389.0_dp,  7062834065.0_dp, -512359803.6_dp, &
         10879881.29_dp,   -86327.92757_dp,  228.4622733_dp, 0.0_dp]
    real(dp), parameter :: rq1(7) = [ &
         40076544269.0_dp, 745249964.8_dp,   7189466.438_dp, &
         47447.26470_dp,   226.1030244_dp,   1.0_dp,         0.0_dp]
    real(dp), parameter :: rp2(5) = [ &
         1.0_dp,           -0.1098628627e-2_dp,  0.2734510407e-4_dp, &
        -0.2073370639e-5_dp, 0.2093887211e-6_dp]
    real(dp), parameter :: rq2(5) = [ &
        -0.1562499995e-1_dp, 0.1430488765e-3_dp, -0.6911147651e-5_dp, &
         0.7621095161e-6_dp,-0.9349451520e-7_dp]
 
    if (x < 8.0_dp) then
      z   = x * x
      res = (rp1(1)+z*(rp1(2)+z*(rp1(3)+z*(rp1(4)+z*(rp1(5)+z*rp1(6)))))) / &
            (rq1(1)+z*(rq1(2)+z*(rq1(3)+z*(rq1(4)+z*(rq1(5)+z*rq1(6)))))) &
            + (2.0_dp/pi) * besj0(x) * log(x)
    else
      z  = 8.0_dp / x
      xx = x - 0.785398164_dp
      res = sqrt(0.636619772_dp/x) * &
            ( sin(xx)*(rp2(1)+z*(rp2(2)+z*(rp2(3)+z*(rp2(4)+z*rp2(5))))) &
            + z*cos(xx)*(rq2(1)+z*(rq2(2)+z*(rq2(3)+z*(rq2(4)+z*rq2(5))))) )
    end if
  end function besy0
 
 
  !---------------------------------------------------------------------------!
  ! besj0(x) : Bessel function of the first kind, order zero
  ! (needed internally by besy0 for the log term)
  !---------------------------------------------------------------------------!
  function besj0(x) result(res)
    real(dp), intent(in) :: x
    real(dp)             :: res
    real(dp)             :: ax, z, xx
 
    real(dp), parameter :: rp1(7) = [ &
        57568490574.0_dp, -13362590354.0_dp, 651619640.7_dp, &
        -11214424.18_dp,   77392.33017_dp,  -184.9052456_dp, 1.0_dp]
    real(dp), parameter :: rq1(7) = [ &
        57568490411.0_dp,  1029532985.0_dp,  9494680.718_dp, &
        59272.64853_dp,    267.8532712_dp,   1.0_dp,         0.0_dp]
    real(dp), parameter :: rp2(5) = [ &
         1.0_dp,           -0.1098628627e-2_dp,  0.2734510407e-4_dp, &
        -0.2073370639e-5_dp, 0.2093887211e-6_dp]
    real(dp), parameter :: rq2(5) = [ &
        -0.1562499995e-1_dp, 0.1430488765e-3_dp, -0.6911147651e-5_dp, &
         0.7621095161e-6_dp,-0.9349451520e-7_dp]
 
    ax = abs(x)
    if (ax < 8.0_dp) then
      z   = x * x
      res = (rp1(1)+z*(rp1(2)+z*(rp1(3)+z*(rp1(4)+z*(rp1(5)+z*(rp1(6)+z*rp1(7))))))) / &
            (rq1(1)+z*(rq1(2)+z*(rq1(3)+z*(rq1(4)+z*(rq1(5)+z*(rq1(6)+z*rq1(7)))))))
    else
      z  = 8.0_dp / ax
      xx = ax - 0.785398164_dp
      res = sqrt(0.636619772_dp/ax) * &
            ( cos(xx)*(rp2(1)+z*(rp2(2)+z*(rp2(3)+z*(rp2(4)+z*rp2(5))))) &
            - z*sin(xx)*(rq2(1)+z*(rq2(2)+z*(rq2(3)+z*(rq2(4)+z*rq2(5))))) )
    end if
  end function besj0
 
 
  !---------------------------------------------------------------------------!
  ! H0(x) : Struve function of order zero
  !
  ! Uses the power series for small x:
  !   H0(x) = (2/pi) * SUM_{m=0}^{inf} (-1)^m * x^(2m+1) / [(2m+1) * (1*3*5...*(2m+1))^... ]
  !
  ! More precisely:
  !   H0(x) = (2/pi) * [ x/1^2 - x^3/(1^2*3^2) + x^5/(1^2*3^2*5^2) - ... ]
  !
  ! For large x, use the asymptotic relation:
  !   H0(x) = Y0(x) + (2/pi) * SUM_{m=0}^{M} (2m-1)!! / x^(2m+1) * (-1)^m ... 
  ! approximated as:
  !   H0(x) ~ Y0(x) + (2/pi)*[1/x + 1*1/(x^3) + 1*1*9/(x^5) + ...]
  !---------------------------------------------------------------------------!
  function struveh0(x) result(res)
    real(dp), intent(in) :: x
    real(dp)             :: res
    real(dp)             :: term, xsq, sign
    integer              :: m
 
    if (x <= 20.0_dp) then
      ! Power series: H0(x) = (2/pi)*sum_{m=0}^inf (-1)^m * x^(2m+1) / [(2m+1)!!]^2
      ! term_{m} = (-1)^m * x^(2m+1) / [(1*3*5*...*(2m+1))^2]
      xsq  = x * x
      term = x                   ! m=0 term: x / 1^2
      res  = term
      sign = -1.0_dp
      do m = 1, 100
        term = term * xsq / real((2*m-1)*(2*m+1), dp)**2 * real((2*m-1),dp)**2
        ! Recurrence: term_m = term_{m-1} * x^2 / (2m+1)^2
        term = -term * xsq / real((2*m+1), dp)**2
        res  = res + term
        if (abs(term) < 1.0e-15_dp * abs(res)) exit
      end do
      res = res * (2.0_dp / pi)
    else
      ! Asymptotic: H0(x) ~ Y0(x) + (2/pi)*[1 - 1/x^2 + 9/x^4 - ...] / x
      ! Use: H0(x) = Y0(x) + (2/pi) * sum_{m=0}^{M} (-1)^m*(2m-1)!!^2 / x^(2m+1)
      ! with (−1)!! = 1 by convention
      term = 1.0_dp / x
      res  = term
      do m = 1, 20
        term = -term * real((2*m-1)*(2*m-1), dp) / (x*x)
        res  = res + term
        if (abs(term) > abs(res)) exit   ! series diverging, stop
      end do
      res = besy0(x) + (2.0_dp/pi) * res
    end if
  end function struveh0
 
end module special_functions_mod
 
 
!=============================================================================!
module quadrature_mod
  use precision_mod
  implicit none
 
contains
 
  !---------------------------------------------------------------------------!
  ! Gauss-Legendre nodes & weights on [-1,1]
  !---------------------------------------------------------------------------!
  subroutine gauss_legendre(n, nodes, weights)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: nodes(n), weights(n)
 
    integer  :: i, j, m
    real(dp) :: gl_p1, gl_p2, gl_p3, pp, z, z1
    real(dp), parameter :: tol = 3.0e-14_dp
 
    m = (n + 1) / 2
    do i = 1, m
      z = cos(pi * (real(i,dp) - 0.25_dp) / (real(n,dp) + 0.5_dp))
      do
        gl_p1 = 1.0_dp;  gl_p2 = 0.0_dp
        do j = 1, n
          gl_p3 = gl_p2;  gl_p2 = gl_p1
          gl_p1 = ((2.0_dp*j-1.0_dp)*z*gl_p2 - (real(j,dp)-1.0_dp)*gl_p3) / real(j,dp)
        end do
        pp = real(n,dp) * (z*gl_p1 - gl_p2) / (z*z - 1.0_dp)
        z1 = z;  z = z1 - gl_p1/pp
        if (abs(z - z1) <= tol) exit
      end do
      nodes(i)       = -z
      nodes(n+1-i)   =  z
      weights(i)     =  2.0_dp / ((1.0_dp - z*z)*pp*pp)
      weights(n+1-i) =  weights(i)
    end do
  end subroutine gauss_legendre
 
 
  !---------------------------------------------------------------------------!
  ! Integrate  e^t * (X^2 + t^2)^(-1/2)  over [0, Y]
  ! using Gauss-Legendre quadrature
  !---------------------------------------------------------------------------!
  function integrate_kernel(X, Y, n_gl) result(Ival)
    real(dp), intent(in) :: X, Y
    integer,  intent(in) :: n_gl
    real(dp)             :: Ival
 
    real(dp), allocatable :: nodes(:), weights(:)
    real(dp) :: mid, half, t
    integer  :: i
 
    allocate(nodes(n_gl), weights(n_gl))
    call gauss_legendre(n_gl, nodes, weights)
 
    ! Map [-1,1] -> [0, Y]
    mid  = 0.5_dp * Y
    half = 0.5_dp * Y
 
    Ival = 0.0_dp
    do i = 1, n_gl
      t    = mid + half * nodes(i)
      Ival = Ival + weights(i) * exp(t) / sqrt(X*X + t*t)
    end do
    Ival = Ival * half
 
    deallocate(nodes, weights)
  end function integrate_kernel
 
end module quadrature_mod
 
 
!=============================================================================!
program fxy_main
  use precision_mod
  use special_functions_mod
  use quadrature_mod
  implicit none
 
  integer,  parameter :: N_GL = 64        ! GL points (high accuracy, no singularity)
 
  ! --- Grid matching the reference table ---
  integer,  parameter :: NY = 8, NX = 9
  real(dp), parameter :: Yvals(NY) = [0.1_dp, 0.2_dp, 0.5_dp, 1.0_dp, &
                                       2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]
  real(dp), parameter :: Xvals(NX) = [0.10_dp, 0.20_dp, 0.50_dp, 1.00_dp, &
                                       2.00_dp, 5.00_dp, 10.00_dp, 20.00_dp, 50.00_dp]
 
  real(dp) :: X, Y, Ikernel, F_val
  integer  :: ix, iy, iunit
 
  open(newunit=iunit, file='fxy_results.csv', status='replace', action='write')
  write(iunit,'(A)', advance='no') 'Y'
  do ix = 1, NX
    write(iunit,'(A,F6.2)', advance='no') ',X=', Xvals(ix)
  end do
  write(iunit,*)
 
  write(*,'(A)') '======================================================================'
  write(*,'(A)') '  F(X,Y) = -2*exp(-Y)*INT_0^Y e^t*(X^2+t^2)^(-1/2) dt'
  write(*,'(A)') '           - pi*exp(-Y)*[Y0(X) + H0(X)]'
  write(*,'(A)') '======================================================================'
 
  do iy = 1, NY
    Y = Yvals(iy)
 
    write(iunit,'(F6.2)', advance='no') Y
    write(*,'(A,F6.2)') '  Y = ', Y
 
    do ix = 1, NX
      X = Xvals(ix)
 
      ! --- Term 1: -2*exp(-Y) * INT_0^Y e^t*(X^2+t^2)^(-1/2) dt ---
      Ikernel = integrate_kernel(X, Y, N_GL)
 
      ! --- Term 2: -pi*exp(-Y)*[Y0(X) + H0(X)] ---
      F_val = -2.0_dp * exp(-Y) * Ikernel &
              - pi * exp(-Y) * (besy0(X) + struveh0(X))
 
      write(*,'(A,F6.2,A,ES22.15)') '    X = ', X, '   F = ', F_val
      write(iunit,'(A,ES22.15)', advance='no') ',', F_val
    end do
 
    write(*,*)
    write(iunit,*)
  end do
 
  close(iunit)
  write(*,'(A)') 'Results written to fxy_results.csv'
 
end program fxy_main

