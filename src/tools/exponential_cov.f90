!> \file exponential_cov.f90
!! Set of tools to build exponential covariance matrices.
!! By now, this toolbox support only the uniformly spaced grid: fixed discretization step.
!!  \author Innocent Souopgui
!!  \date 12/2013
!!  \version 1.0
!<


!> Procedures to build exponential covariance matrices
!!
!<
module exponential_cov
  use general_constant
  use general_tools
  use debug_tools

implicit none

  private
  public exp_cov_param, se_cov_matrix, load_ecp

  !> Derived type for 1D exponential covariance matrix parameters
  type exp_cov_param
    integer:: n !< Number of grid point
    integer:: k !< Number of significant eigenvalues/singular values to be kept
    real(dp):: d !< discretization step
    real(dp):: s !< smoothing parameter of the covariance function (variance at a single point)
    real(dp):: l !< decorrelation length
  end type exp_cov_param

  !> Computes a square exponential covariance matrix
  interface se_cov_matrix
    module procedure se_cov_matrix_scalar
    module procedure se_cov_matrix_type
  end interface se_cov_matrix

contains

  !> \brief Loads the parameters of an exponential covariance matrix from namelist file
  !! \param[in] ecp data structure containing the exponential covariance matrix parameters
  !! \param[in] fName name of the namelist file
  !! \param[in] dName name associated with physical space dimension where the covariance matrix will be applied.
  !! \param[in] n size of the covariance matrix (number of columns/rows)
  !! \param[in] delta discretization step
  !!
  !! \details The namelist is supposed to holds
  !! - the dimension name to make the difference with other covariance parameters in the namelist file
  !! - the number of significant singular values
  !! - the smoothing parameter
  !! - and the decorrelation length
  !! By default, dName picks values in ["X", "Y", "Z"], but can be anything meaningful to the program writer
  !<
  subroutine load_ecp(ecp, fName, dName, n, delta)
    INTEGER, PARAMETER:: fid = 68
    type(exp_cov_param), intent(out) :: ecp
    CHARACTER(len=*), INTENT(IN) :: fName, dName
    INTEGER, INTENT(IN) :: n
    REAL(cp), INTENT(IN) :: delta
    !local var
    CHARACTER(len=ip_snl) :: dimName, ala_dimName
    INTEGER :: nsv
    REAL(cp):: smoothing, decorrelation
    namelist/NAM_EXPONENTIAL_COV_PARAM/&
      dimName,&
      nsv,&
      smoothing,&
      decorrelation

    open(fid,file=fName,form='FORMATTED',status='OLD')
    dimName = "{#}@"
    ala_dimName = trim(dName)
    call uppercase(ala_dimName)
    call debug(ala_dimName, "In load_ecp, looking for dimName = ", tag=dALLWAYS)
    do while(dimName.ne.ala_dimName)
      read(fid, NAM_EXPONENTIAL_COV_PARAM)!reading the block
      call uppercase(dimName)
      call debug(dimName, "  found dimName = ", tag=dALLWAYS)
    end do
    close(fid)

    ecp = exp_cov_param(n, nsv, delta, smoothing, decorrelation)
  end subroutine load_ecp

  !> Computes the square exponential covariance matrix for a 1D uniform grid, using exp_cov_param derived type
  !! @param [in, out] C square matrix of covariance
  !! @param [in] xcp covariance parameters in the x dimension
  !!
  !! @see se_cov_matrix_scalar
  !>
  subroutine se_cov_matrix_type(C, xcp)
    real(dp), dimension(:,:), intent(in out) :: C
    type(exp_cov_param), intent(in) :: xcp

    call se_cov_matrix_scalar(C, xcp%d, xcp%s, xcp%l)
  end subroutine se_cov_matrix_type

  !> Computes the square exponential covariance matrix for a 1D uniform grid, using scalar parameters
  !! @param [in, out] C square matrix of covariance
  !! @param [in] dx discretization step
  !! @param [in] s smoothing parameter of the covariance function
  !! @param [in] lx correlation length
  !!
  !! The covariance betwen the grid points i and j is given by
  !! \f$ C(i,j) = s^2\exp(\frac{-(i\Delta x-j\Delta x)^2}{2l_x^2}) \f$
  !>
  subroutine se_cov_matrix_scalar(C, dx, s, lx)
    real(dp), dimension(:,:), intent(in out) :: C
    real(dp), intent(in) :: dx, s, lx
    !local variables
    real(dp) :: x, s2, lx2
    integer :: ii, jj

    s2 = s*s
    lx2 = lx*lx

    !compute the upper triangular part
    do jj=2, size(C,2)
      do ii=1, jj-1
        x = (ii-jj)*dx
        C(ii, jj) = s2*exp(-x*x/(2*lx2))
      end do
    end do

    !The diagonal element are all identical
    do jj=1,size(C,2)
      C(jj, jj) = s2
    end do

    !computes the lower triangular part by symmetry
    do jj = 1, size(C,2)
      do ii = jj+1, size(C,1)
        C(ii, jj) = C(jj, ii)
      end do
    end do

    !I strongly believe that the next (commented line) was not necessary
    !C = dx*C
  end subroutine se_cov_matrix_scalar

  !> Covariance function using exponential of the square
  !! @param [in] x1, x2 points between which the covariance is evaluated
  !! @param [in] s smoothing parameter of the covariance function
  !! @param [in] lx decorrelation length
  !!
  !<
  function se_cov_function(x1, x2, s, lx)result(f)
    real(dp), intent(in) :: x1, x2, s, lx
    !local variables
    real(dp) :: f, x

    x = x1-x2
    f = s*s*exp(-x*x/(2*lx*lx))
  end function se_cov_function

  !> Covariance function using exponential of absolute value
  !! @param [in] x1, x2 points between which the covariance is evaluated
  !! @param [in] s smoothing parameter of the covariance function; this is the varaince
  !! @param [in] lx decorrelation length
  !!
  !<
  function ae_cov_function(x1, x2, s, lx)result(f)
    real(dp), intent(in) :: x1, x2, s, lx
    real(dp) :: f

    f = s*s*exp(-abs(x1-x2)/lx);
  end function ae_cov_function
end module exponential_cov