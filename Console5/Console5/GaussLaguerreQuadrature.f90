 module GAUSS_LAGUREE_MODULE
    
    use GAUSS_BASIC_MODULE
    implicit none
   
    type, EXTENDS( GAUSS_BASIC_TYPE )  :: GAUSS_LAGUERRE_TYPE
  
      
    contains   
       procedure   ::   CALCUATE_SEEDS_WEIGHTS  => CALCUATE_SW_LAGUERRE
 
 
    procedure    ::  SORT_PARIS
    procedure    ::  CALCULATE_TRIDIAG

 
   
    end type GAUSS_LAGUERRE_TYPE
    
    
 
    
    
    contains
    
    subroutine CALCUATE_SW_LAGUERRE( THIS )
    implicit none
    
    class( GAUSS_LAGUERRE_TYPE), intent( inout ) :: THIS
    
    !-------------------LOCAL VARIABLES------------------!
    integer   ::   ALLOCATE_ARRAY_FLAG
    integer   ::   II, INFOR_FLAG
    
    real(8)  ::  mu0,TEMP_A
    
    real(8), dimension(:), allocatable    ::  SEEDS_POINT_TEMP, EE_TEMP   
    real(8), dimension(:,:), allocatable  ::  ZZ_TEMP
     
    allocate(   SEEDS_POINT_TEMP(THIS%NUMBER_POINTS), EE_TEMP(THIS%NUMBER_POINTS), &
                 ZZ_TEMP(THIS%NUMBER_POINTS,THIS%NUMBER_POINTS), stat = ALLOCATE_ARRAY_FLAG )
    if( ALLOCATE_ARRAY_FLAG  /= 0)then
             write(*,*)  'ERROR IN ALLCOATE THE ARRAY INSIDE GAUSS LAGURRE MODULE'
        PAUSE   
    end if
    

    ! Build Jacobi matrix (symmetric tridiagonal) for generalized Laguerre weight x^a e^{-x}.
    ! Diagonal: d_i = 2i - 1 + a
    ! Off-diagonal: e_i = sqrt(i*(i+a)), with e_n unused (set to 0)
    
    TEMP_A = 0.0D0
    
    do II  = 1, THIS%NUMBER_POINTS
        
      SEEDS_POINT_TEMP( II ) = 2d0*II - 1.0D0 + TEMP_A
      
    end do
    
    
    do II = 1, THIS%NUMBER_POINTS-1
        
      EE_TEMP( II ) = sqrt( dble(II) * (dble(II) + TEMP_A) )
      
    end do
    
    
    EE_TEMP(THIS%NUMBER_POINTS) = 0.0d0

    ! z initialized to identity; will contain eigenvectors
    ZZ_TEMP(:,:) = 0.0D0
    
    do II = 1, THIS%NUMBER_POINTS
        
      ZZ_TEMP(II,II) = 1.0D0

    end do

    call THIS%CALCULATE_TRIDIAG(SEEDS_POINT_TEMP, EE_TEMP, THIS%NUMBER_POINTS, ZZ_TEMP, INFOR_FLAG)
    
    if (INFOR_FLAG /= 0) then
        
      error stop "gauss_laguerre: eigen-solver failed to converge."
      
    end if

    ! Nodes are eigenvalues
    THIS%SEEDS_POINT(:) = SEEDS_POINT_TEMP(:)

    ! μ0 = ∫_0^∞ x^a e^{-x} dx = Γ(a+1)
    mu0 = gamma(TEMP_A + 1.0D0)

    ! Weights: μ0 * (first component of eigenvector)^2
    do II = 1, THIS%NUMBER_POINTS
        
      THIS%WEIGHTS_POINT(  II ) = mu0 * (ZZ_TEMP(1,II) * ZZ_TEMP(1,II))
      
    end do

    call THIS%SORT_PARIS 

    deallocate(SEEDS_POINT_TEMP , EE_TEMP ,  ZZ_TEMP , stat = ALLOCATE_ARRAY_FLAG )
    
    if( ALLOCATE_ARRAY_FLAG  /= 0)then
        
        write(*,*)  'ERROR IN DEALLCOATE THE ARRAY INSIDE GAUSS LAGURRE MODULE'
        PAUSE
    end if
    
    THIS%SEEDS_FORMING_FLAG  = .TRUE. 
    
    
    end subroutine CALCUATE_SW_LAGUERRE
    
    
    subroutine SORT_PARIS( THIS )
     class( GAUSS_LAGUERRE_TYPE), intent( inout ) :: THIS
     
    integer :: II, JJ
    real(8) :: SEED_TEMP, WEIGHT_TEMP
    
     
    
    do II = 1, THIS%NUMBER_POINTS - 1
        
      do JJ = II + 1, THIS%NUMBER_POINTS
          
        if (   THIS%SEEDS_POINT( JJ ) < THIS%SEEDS_POINT( II )  ) then
            
          SEED_TEMP = THIS%SEEDS_POINT( II )
          THIS%SEEDS_POINT( II ) = THIS%SEEDS_POINT( JJ ) 
          THIS%SEEDS_POINT( JJ ) = SEED_TEMP
          
          
          WEIGHT_TEMP = THIS%WEIGHTS_POINT(  II )
          THIS%WEIGHTS_POINT( II ) = THIS%WEIGHTS_POINT( JJ )
          THIS%WEIGHTS_POINT( JJ ) = WEIGHT_TEMP
        end if
        
      end do
      
    end do
     
     
    
    
    end  subroutine SORT_PARIS
    
    subroutine CALCULATE_TRIDIAG(THIS, d, e, n, z, info)
    
    implicit none
    class( GAUSS_LAGUERRE_TYPE), intent( inout ) :: THIS
    
    
    integer, intent(in)    :: n
    real(8), intent(inout) :: d(n), e(n), z(n,n)
    integer, intent(out)   :: info

    integer :: l, m, i, k, iter
    real(8) :: s, r, p, g, f, dd, c, b

    info = 0

    ! Shift e down (NR convention), but we already have e(1..n-1), e(n)=0.
    ! We'll keep as-is and use e(m) as subdiagonal.

    do l = 1, n
      iter = 0
10    continue
      ! Find small subdiagonal element
      do m = l, n-1
        dd = abs(d(m)) + abs(d(m+1))
        if (abs(e(m)) <= 1d-16*dd) exit
      end do

      if (m /= l) then
        iter = iter + 1
        if (iter > 60) then
          info = 1
          return
        end if

        g = (d(l+1) - d(l)) / (2d0*e(l))
        r = sqrt(g*g + 1d0)
        g = d(m) - d(l) + e(l) / (g + sign(r, g))
        s = 1d0
        c = 1d0
        p = 0d0

        do i = m-1, l, -1
          f = s*e(i)
          b = c*e(i)
          if (abs(f) >= abs(g)) then
            c = g/f
            r = sqrt(c*c + 1d0)
            e(i+1) = f*r
            s = 1d0/r
            c = c*s
          else
            s = f/g
            r = sqrt(s*s + 1d0)
            e(i+1) = g*r
            c = 1d0/r
            s = s*c
          end if
          g = d(i+1) - p
          r = (d(i) - g)*s + 2d0*c*b
          p = s*r
          d(i+1) = g + p
          g = c*r - b

          ! Apply rotation to eigenvector matrix z
          do k = 1, n
            f       = z(k,i+1)
            z(k,i+1)= s*z(k,i) + c*f
            z(k,i)  = c*z(k,i) - s*f
          end do
        end do

        d(l) = d(l) - p
        e(l) = g
        e(m) = 0d0
        goto 10
      end if
    end do

  end subroutine CALCULATE_TRIDIAG
    
    end module GAUSS_LAGUREE_MODULE