!=============================================================================!
!  green_transform_mod.f90
!
!  Transforms F(X,Y) and its derivatives (as returned by integrate_kernel,
!  in terms of the reduced radial/vertical coordinates X, Y) into the
!  physical Green's function G and its derivatives w.r.t. the real
!  Cartesian coordinates (x,y,z) and image-point coordinates (zeta,eta,xi).
!
!  Geometry (from problem definition):
!
!     X = sqrt[ (xbar-zetabar)^2 + (ybar-etabar)^2 ]
!     Y = -( zbar + xibar )
!
!     xbar = k0*x,  ybar = k0*y,  zbar = k0*z
!     zetabar = k0*zeta,  etabar = k0*eta,  xibar = k0*xi
!
!     G = k0 * F(X,Y)
!
!  First derivatives (as given):
!
!     dG/dx    = k0^2 * xhat * dF/dX      ,  dG/dzeta = -dG/dx
!     dG/dy    = k0^2 * yhat * dF/dX      ,  dG/deta  = -dG/dy
!     dG/dz    = -k0^2 * dF/dY            ,  dG/dxi   =  dG/dz
!
!  where xhat = (xbar-zetabar)/X ,  yhat = (ybar-etabar)/X   (direction
!  cosines of the horizontal offset vector).
!
!  Second derivatives (NOT shown in the original figure -- derived here by
!  applying the same chain rule a second time; included because your
!  integrate_kernel already returns F_XX, F_XY, F_YY. Delete this block if
!  you only need the first-derivative transform).
!
!     Let u = xbar-zetabar, v = ybar-etabar, so X = sqrt(u^2+v^2).
!
!     d2G/dx2   = k0^3 * [ (v^2/X^3)*F_X + (u^2/X^2)*F_XX ]
!     d2G/dy2   = k0^3 * [ (u^2/X^3)*F_X + (v^2/X^2)*F_XX ]
!     d2G/dxdy  = k0^3 * (u*v/X^2) * [ F_XX - F_X/X ]
!     d2G/dz2   = k0^3 * F_YY
!     d2G/dxdz  = -k0^3 * (u/X) * F_XY
!     d2G/dydz  = -k0^3 * (v/X) * F_XY
!
!     (Mixed derivatives w.r.t. the image coordinates zeta,eta,xi follow
!      from the same antisymmetry used for the first derivatives, e.g.
!      d2G/dzeta2 = d2G/dx2, d2G/dxdzeta = -d2G/dx2, etc. -- add as needed.)
!=============================================================================!

module GREEN_TRANSFORM_MOD

  use precision_mod
  implicit none

  private
  public :: TRANSFORM_GREEN_DERIVATIVES

  real(DP), parameter :: SAFE_MIN_RADIUS = 1.0E-14_DP

contains

  !---------------------------------------------------------------------!
  ! TRANSFORM_GREEN_DERIVATIVES
  !
  ! Inputs:
  !   X, Y, Z, ZETA, ETA, XI  - physical / image coordinates
  !   K0                      - wavenumber (scaling factor)
  !   F, FX, FY, FXX, FXY, FYY  - F(X,Y) and its X/Y derivatives, i.e.
  !                               KERNEL_OUT(1:6) from integrate_kernel
  !
  ! Outputs:
  !   G                      - k0*F
  !   GX, GY, GZ             - dG/dx, dG/dy, dG/dz
  !   GZETA, GETA, GXI       - dG/dzeta, dG/deta, dG/dxi
  !   GXX, GYY, GXY,
  !   GZZ, GXZ, GYZ          - optional second derivatives (real space)
  !---------------------------------------------------------------------!
  subroutine TRANSFORM_GREEN_DERIVATIVES( X, Y, Z, ZETA, ETA, XI, K0, &
                                           F, FX, FY, FXX, FXY, FYY,  &
                                           G, GX, GY, GZ, GZETA, GETA, GXI, &
                                           GXX, GYY, GXY, GZZ, GXZ, GYZ )

    real(DP), intent(in)  :: X, Y, Z, ZETA, ETA, XI, K0
    real(DP), intent(in)  :: F, FX, FY, FXX, FXY, FYY

    real(DP), intent(out) :: G, GX, GY, GZ, GZETA, GETA, GXI
    real(DP), intent(out), optional :: GXX, GYY, GXY, GZZ, GXZ, GYZ

    real(DP) :: XBAR, YBAR, ZBAR, ZETABAR, ETABAR, XIBAR
    real(DP) :: U, V, XR, XHAT, YHAT
    real(DP) :: K0SQ, K0CU

    !--- scaled coordinates -------------------------------------------!
    XBAR    = K0 * X
    YBAR    = K0 * Y
    ZBAR    = K0 * Z
    ZETABAR = K0 * ZETA
    ETABAR  = K0 * ETA
    XIBAR   = K0 * XI

    U  = XBAR - ZETABAR
    V  = YBAR - ETABAR
    XR = sqrt( U**2 + V**2 )

    !--- guard against the singular point X -> 0 -----------------------!
    if ( XR < SAFE_MIN_RADIUS ) XR = SAFE_MIN_RADIUS

    XHAT = U / XR
    YHAT = V / XR

    K0SQ = K0 * K0
    K0CU = K0SQ * K0

    !--- function value --------------------------------------------------!
    G = K0 * F

    !--- first derivatives -------------------------------------------!
    GX = K0SQ * XHAT * FX
    GY = K0SQ * YHAT * FX
    GZ = -K0SQ * FY

    GZETA = -GX
    GETA  = -GY
    GXI   =  GZ

    !--- optional second derivatives (real-space Hessian) --------------!
    if ( present(GXX) ) then
        GXX = K0CU * ( (V**2 / XR**3) * FX + (U**2 / XR**2) * FXX )
    end if
    if ( present(GYY) ) then
        GYY = K0CU * ( (U**2 / XR**3) * FX + (V**2 / XR**2) * FXX )
    end if
    if ( present(GXY) ) then
        GXY = K0CU * ( U * V / XR**2 ) * ( FXX - FX / XR )
    end if
    if ( present(GZZ) ) then
        GZZ = K0CU * FYY
    end if
    if ( present(GXZ) ) then
        GXZ = -K0CU * ( U / XR ) * FXY
    end if
    if ( present(GYZ) ) then
        GYZ = -K0CU * ( V / XR ) * FXY
    end if

  end subroutine TRANSFORM_GREEN_DERIVATIVES

end module GREEN_TRANSFORM_MOD
