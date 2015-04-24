!> @file evd_covariance.f90
!! Defines a general module for covariance matrices.
!<

!> @brief General module for covariance matrix
!!
!! @author Innocent Souopgui
!!
!! \details This module implements a set of subroutines for the covariance
!!   matrix for one and two dimensional problems
!!  Using the positive (semi) definiteness of the covariance matrices,
!!  the covariance matrix B is diagonalizable and its singular value
!!  decomposition is also its diagonalization in a basis of eigenvectors.
!! We consider the case where B is approximated using only the leading
!! eigenvectors, the full case being similar.
!! The covariance matrix is supposed to be given under a factorized forms
!!   \f[ B \approx UDU^T \f] where the colums of U are the leading eigenvectors
!!   of B and D is the diagonal matrix of the leading eigenvalues of B.
!! For the 2D, we consider the covariance function that factorises as a
!! product of covariance functions in each spatial dimension; so that the
!! covariance matrix B is the kronecker product of the covariance  matrices
!! in the single dimensions \f[ B = Kron(Bx,By) \f]. Bx and By are given under
!! their factorized forms \f[ B = Kron(Ux,Uy)*Kron(Dx,Dy)*Kron(Ux,Uy)^T \f],
!! where Dx is a diagonal matrix with leading eigenvalues of Bx on the diagonal
!! and Dy is a diagonal matrix with eigenvalues of By on the diagonal
!! Bx = Ux*Dx*Ux^T,  By = Uy*Dy*Uy^T.
!!The 2D covariance matrix applies to a vector that correspond to a discretization
!!  of a 2d function (function of two variables) on the same grid as the control
!!  variable. The chronecker product factorization depends on the ordering of
!!  the dimensions in the discretization as well as the vectorization of the 2D
!!  discretized function.
!!Suppose that the covariance function factorizes such that the covariances on
!!  the grid are combinations of the covariances in each direction,
!!  cov(p1, p2) = cov_x(x1,x2)cov_y(y1,y2) where p1 is the point of coordinates
!!  (x1,y1) and p2 (x2,y2). In this case B*f = Kron(B_x,B_y)*f where B_x and B_y
!!  are respectively the covariances in the direction x, resp. y.
!! According to the ordering of the discretization points and the vectorization
!!  order, that formula can change:
!! Let F be the matrix of a discretized function, and f its vectorized
!! representation.
!! 1) if the first dimension of F corresponds to the x direction (x-leading) then:
!!   - g=B*f = Kron(B_x,B_y)*f, if f is a row major vectorization
!!     (C-prefered ordering);
!!   - g=B*f = Kron(B_y,B_x)*f, if f is a column major vectorization
!!       (Fortran-prefered ordering).
!! 2) if the first dimension of F corresponds to the y direction (y-leading) then:
!!   - g=B*f = Kron(B_y,B_x)*f with a row major vectorization
!!      (C-prefered ordering) of F
!!   - g=B*f = Kron(B_x,B_y)*f, with a column major vectorization
!!      (Fortran-prefered ordering) of F.
!! All the 2D routines in this module are based on the formula
!! g=B*f = Kron(B_x,B_y)*f. It is assumed either :
!! - x-leading and row major vectorization or
!! - y-leading and column major vectorization
!!
!<
MODULE evd_covariance
  USE general_constant
  USE general_tools, only: readInfo, readMatrix, writeMatrix, fileExist
  USE debug_tools
  use kronecker
  use exponential_cov
  use evd_tools
IMPLICIT NONE

  private
  public evd_cov_matrix&
    , load_evd_covariance, save_evd_covariance&
    , apply_evd_cov, apply_evd_cov_inv&
    , apply_evd_square_cov, apply_evd_square_cov_inv&
    , apply_evd_square_covb&
    , get_evd_cov_nsv, make_evd_cov, get_evd_cov_size&
    , get_evd_cov_nRow
  !, cov_mv, cov_inv_mv, cov_square_mv, cov_square_mvb, cov_square_inv_mv&


  !> \brief default number of factors in the Kronecker factorization
  !!  of the covariance matrix, used for non initialized covariance matrix
  !<
  INTEGER, PARAMETER :: default_nFact=-1

    !> Derived type for singular values decomposed covariance matrices
    !! It is assumed that the covariance matrix is provided as its Singular
    !! Value Decomposition: B = UDUt. Furthermore, it can be (Kronecker)
    !! factorized along the physical dimensions of the problem:
    !! B = Bx kron By kron Bz
    type evd_cov_matrix

        !> \brief number of terms in the Kronecker factorization of
        !!  the covariance matrix.
        !!Conceptually, it corresponds to the number of physical dimensions
        !! of the problem if the covariance is factorized.
        !!Logically, it can be any Kronecker factorization.
        !<
        INTEGER, private :: im_nFac = default_nFact

        !> \brief matrix(ces) of singular vector of the covariance matrix
        !! \details If the covariance matrix B is not factorized, only Ux is initialized
        !! IF B is factorized as B = Bx kron By (two terms), Ux and Uy are provided
        !! A maximum of 3 terms in the factorization is supported.
        !<
        REAL(cp), private, ALLOCATABLE, DIMENSION(:,:):: Ux, Uy, Uz

        !> vector(s) of eigenvalues of B
        REAL(cp), private, ALLOCATABLE, DIMENSION(:) :: Dx, Dy, Dz

    end type evd_cov_matrix

!private procedures interface
  !> Covariance matrix vector product
  interface cov_mv
    module procedure cov_mv_default
    module procedure cov_mv_2fac
    module procedure cov_mv_3fac
  end interface cov_mv

  !> Inverse covariance matrix vector product
  interface cov_inv_mv
    module procedure cov_inv_mv_default
    module procedure cov_inv_mv_2fac
    module procedure cov_inv_mv_3fac
  end interface cov_inv_mv

  !> square root covariance matrix vector product
  !! Those subroutines actually apply the operator L derived from
  !! the eigenvalue decomposition of the covariance matrix.
  !! let B = UDU^t be the eigenvalue decomposition
  !! Those subroutines apply UD^(1/2)
  !! To apply the square root that should be UD^(1/2)U^T
  !! Those subroutines are used to generate multivariate random variable
  !! with covariance B from the standard normal distribution.
  !! Since any decomposition of the form LL^t can be used.
  !! This should be corrected in the future to apply the square root
  !! That can also be used for the multivariate random numbers generation.
  !! It has an overhead computation with the application of U^t.
  !<
  interface cov_square_mv
    module procedure cov_square_mv_default
    !module procedure cov_square_mv_1d
    module procedure cov_square_mv_2fac
    module procedure cov_square_mv_3fac
  end interface cov_square_mv

  !> square root covariance matrix vector product (adjoint)
  interface cov_square_mvb
    module procedure cov_square_mvb_default
    module procedure cov_square_mvb_2fac
    !module procedure cov_square_mv_3fac
  end interface cov_square_mvb

  !> Inverse square root covariance matrix vector product
  interface cov_square_inv_mv
    module procedure cov_square_inv_mv_default
    module procedure cov_square_inv_mv_2fac
    module procedure cov_square_inv_mv_3fac
  end interface cov_square_inv_mv

  !> Checks the parameters for the application of the covariance matrix
  !! or its inverse
  interface cov_mv_check
    module procedure cov_mv_check_default
    module procedure cov_mv_check_2fac
    module procedure cov_mv_check_3fac
  end interface cov_mv_check

  !> Check the parameters for the application of the square root of the
  !! covariance matrix or its inverse
  interface cov_square_mv_check
    module procedure cov_square_mv_check_default
    module procedure cov_square_mv_check_2fac
    module procedure cov_square_mv_check_3fac
  end interface cov_square_mv_check
!end of private procedures interface

contains

!public procedures

  !> \brief returns the number of singular values used to define the
  !!  approximation of the covariance matrix; that is also the size of
  !!  the changed state variable
  !! @param [in] C covariance matrix derived type
  !! @param [in] [nctl] size of the control vector
  !! @details
  !! If a call to load_covariance has not been done before: the function
  !! returns nctl if it is provided, otherwise, the program is stopped.
  !! If a call to load_covariance have been done before, the associated value
  !! is returned
  !<
  FUNCTION get_evd_cov_nsv(C, nctl) RESULT(nsv)
    type(evd_cov_matrix), intent(in) :: C
    INTEGER, OPTIONAL, INTENT(IN) :: nctl
    !local variables
    INTEGER :: nsv

    SELECT CASE(C%im_nFac)
      CASE (1)
        nsv = SIZE(C%Dx)
      CASE (2)
        nsv = SIZE(C%Dx)*SIZE(C%Dy)
      CASE (3)
        nsv = SIZE(C%Dx)*SIZE(C%Dy)*SIZE(C%Dz)
      CASE (default_nFact)!default value
        IF( PRESENT(nctl) )THEN
          nsv = nctl
        ELSE
          CALL stop_program( 'In evd_cov_nsv, covariance not loaded and nctl not provided' )
        END IF
      CASE DEFAULT
        CALL stop_program( 'In evd_cov_nsv, covariance operations are define only for 1D and 2D' )
    END SELECT
  END FUNCTION get_evd_cov_nsv

  !> \brief returns the number of rows of the covariance matrix
  !! @param [in] C covariance matrix derived type
  !!
  !! \details If \a C is not initialized, the program program stops.
  !<
  FUNCTION get_evd_cov_nRow(C) RESULT(nRow)
    type(evd_cov_matrix), intent(in) :: C
    !local variables
    INTEGER :: nRow

    SELECT CASE(C%im_nFac)
      CASE (1)
        nRow = SIZE(C%Ux, 1)
      CASE (2)
        nRow = SIZE(C%Ux, 1)*SIZE(C%Uy, 1)
      CASE (3)
        nRow = SIZE(C%Ux, 1)*SIZE(C%Uy, 1)*SIZE(C%Uz, 1)
      CASE DEFAULT
        CALL stop_program( 'In get_evd_cov_nRow, covariance not initialized' )
    END SELECT
  END FUNCTION get_evd_cov_nRow

  !> \brief Reads a covariance matrix stored under the form of an evd decomposition
  !! @param [in] UFName file name of the matrix of singular vectors
  !! @param [in] DFName file name of the matrix of singular values
  !! @param [in, out] U matrix of singular vectors
  !! @param [in, out] D vector of singular values
  !<
  SUBROUTINE readEVD( UFName, DFName, U, D )
    CHARACTER(len=*), INTENT(IN)  :: UFName, DFName
    REAL(cp), DIMENSION(:,:), INTENT(INOUT) :: U
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: D
    !local variables
    REAL(cp), DIMENSION(:,:), ALLOCATABLE :: T
    INTEGER :: il_n, il_k, il_nRow, il_nCol, ibj

    CALL readEVDSize( UFName, DFName, il_n, il_k, il_nRow, il_nCol )
    ALLOCATE( T(il_nRow, il_nCol) )

    CALL readMatrix( UFName, U )
    CALL readMatrix( DFName, T )

    IF(il_nRow == 1) THEN
      D = T(1,:)
    ELSE IF(il_nCol == 1) THEN
      D = T(:,1)
    ELSE
      DO ibj = 1, il_nCol
        D(ibj) = T( ibj, ibj )
      END DO
    END IF
  END SUBROUTINE readEVD

  !> \brief read the size of the EV decomposition of B: size of each singular
  !!   vector, and number of singular vector
  !! @param [in] UFName file name of the matrix of singular vectors
  !! @param [in] DFName file name of the matrix of singular values
  !! @param [in, out] n size of each singular vector (number of rows/colums of B)
  !! @param [in, out] k number of singular values (or vector)
  !! @param [in, out] nrD (optional) number of rows of D
  !! @param [in, out] ncD (optional) number of colums of D
  !! if the sizes in UFName and DFName do not match, the program terminates
  !<
  SUBROUTINE readEVDSize(UFName, DFName, n, k, nrD, ncD)
    CHARACTER(len=*), INTENT(IN)  :: UFName, DFName
    INTEGER, INTENT(OUT) :: n, k
    INTEGER, OPTIONAL, INTENT(OUT) :: nrD, ncD
    INTEGER :: il_nrD, il_ncD

    CALL readInfo( UFName, n, k )
    CALL readInfo( DFName, il_nrD, il_ncD )
    CALL EVD_check( n, k, il_nrD, il_ncD )
    IF( PRESENT(nrD) ) nrD = il_nrD
    IF( PRESENT(ncD) ) ncD = il_ncD
  END SUBROUTINE readEVDSize

  !> \brief Loads the covariance matrix represented by its singular value
  !!   decomposition
  !! @param [in, out] C covariance matrix derived type
  !! @param [in] dir directory containing the data files
  !!
  !! \details It is assumed that the covariance matrix B is provided as
  !! its Singular Value Decomposition: B = UDU^t. Furthermore, it can
  !! be (Kronecker) factorized (for example along the physical dimensions of
  !! the problem) in which case it is given under the form:
  !! - B = Bx kron By = (UxDxUx^t) kron (UyDyUy^t) or
  !! - B = Bx kron By kron Bz = (UxDxUx^t) kron (UyDyUy^t) kron (UzDzUz^t)
  !! Conventions for the file names. Data files are named after the matrices
  !! names with the extension 'dat'
  !! - 'U.dat' and 'D.dat' for the non factorized case
  !! - 'Ux.dat',  'Uy.dat', 'Uz.dat' and 'Dx.dat', 'Dy.dat', 'Dz.dat' for
  !!   the factorized case
  !! the existence of files is checked to determine the appropriate configuration.
  !! ('U.dat' and 'D.dat') or ('Ux.dat' and 'Dx.dat') should be present,
  !! otherwise the program terminates
  !!
  !<
  SUBROUTINE load_evd_covariance(C, dir)
    type(evd_cov_matrix), intent(in out) :: C
    CHARACTER(len=*), INTENT(IN)  :: dir
    !local variables
    INTEGER :: il_N, il_K
    CHARACTER(len=ip_fnl) :: UFName, DFName

    CALL debug('Entering load_covariance ', tag=dALLWAYS )
    !Checking files
    IF( fileExist(TRIM(dir)//'/U.dat'))THEN
      UFName = TRIM(dir)//'/U.dat'
      DFName = TRIM(dir)//'/D.dat'
    ELSE
      UFName = TRIM(dir)//'/Ux.dat'
      DFName = TRIM(dir)//'/Dx.dat'
    END IF

    CALL debug(100, ' In load_evd_covariance ', tag=dALLWAYS )
    CALL readEVDSize(UFName, DFName, il_N, il_K)
    ALLOCATE( C%Ux(il_N, il_K), C%Dx(il_K) )
    CALL readEVD( UFName, DFName, C%Ux, C%Dx )
    CALL debug(C%Ux, 'Ux = ')
    CALL debug(C%Dx, 'Dx = ')
    C%im_nFac = 1

    CALL debug(200, ' In load_evd_covariance ', tag=dALLWAYS )
    UFName = TRIM(dir)//'/Uy.dat'
    DFName = TRIM(dir)//'/Dy.dat'
    IF( fileExist(uFName) )THEN
      CALL readEVDSize(UFName, DFName, il_N, il_K)
      ALLOCATE( C%Uy(il_N, il_K), C%Dy(il_K) )
      CALL readEVD( UFName, DFName, C%Uy, C%Dy )
      CALL debug(C%Uy, 'Uy = ')
      CALL debug(C%Dy, 'Dy = ')
      C%im_nFac = C%im_nFac + 1
    END IF

    CALL debug('Exiting load_evd_covariance ', tag=dALLWAYS )
  END SUBROUTINE load_evd_covariance


  !> \brief Writes a covariance matrix under the form of an eigenvalues
  !!  decomposition
  !! @param [in] UFName file name of the matrix of singular vectors
  !! @param [in] DFName file name of the vector(matrix) of singular values
  !! @param [in, out] U matrix of singular vectors
  !! @param [in, out] D vector of singular values
  !<
  subroutine writeEVD( UFName, DFName, U, D )
    character(len=*), intent(in)  :: UFName, DFName
    real(cp), dimension(:,:), intent(in) :: U
    real(cp), dimension(:), intent(in) :: D
    !local variables
    real(cp), dimension(size(D),1) :: T

    T(:,1) = D

    call writeMatrix( UFName, U )
    call writeMatrix( DFName, T )
  end subroutine writeEVD

  !> \brief Saves the covariance matrix represented by its eigenvalues
  !!  decomposition
  !! @param [in] C covariance matrix derived type
  !! @param [in] dir directory where the data files will be created
  !! @see load_evd_covariance
  !<
  subroutine save_evd_covariance(C, dir)
    type(evd_cov_matrix), intent(in) :: C
    character(len=*), intent(in)  :: dir
    !local variables
    character(len=ip_fnl) :: UFName, DFName

    call debug('Entering save_evd_covariance ', tag=dALLWAYS )
    !Checking files
    if( C%im_nFac == 1 )then
      UFName = TRIM(dir)//'/U.dat'
      DFName = TRIM(dir)//'/D.dat'
    else if ( C%im_nFac > 1 )then
      UFName = TRIM(dir)//'/Ux.dat'
      DFName = TRIM(dir)//'/Dx.dat'
    else if ( C%im_nFac < 1 )then
      call stop_program('can not save uninitialized covariance matrix')
    end if

    call debug(100, ' In save_evd_covariance ', tag=dALLWAYS )
    call writeEVD( UFName, DFName, C%Ux, C%Dx )

    call debug(200, ' In save_evd_covariance ', tag=dALLWAYS )
    if( C%im_nFac>=2 )then
      UFName = TRIM(dir)//'/Uy.dat'
      DFName = TRIM(dir)//'/Dy.dat'
      call writeEVD( UFName, DFName, C%Uy, C%Dy )
    end if

    call debug(300, ' In save_evd_covariance ', tag=dALLWAYS )
    if( C%im_nFac>=3 )then
      UFName = trim(dir)//'/Uz.dat'
      DFName = trim(dir)//'/Dz.dat'
      call writeEVD( UFName, DFName, C%Uz, C%Dz )
    end if

    call debug('Exiting save_evd_covariance ', tag=dALLWAYS )
  end subroutine save_evd_covariance

  !> Gets the sizes of the factorized arrays that make up the EV decomposed
  !! covariance matrix
  !! @param [in] C covariance matrix derived type
  !! @param [out] Nx number of components [in the x direction]
  !! @param [out] Kx number of singular vectors [in the x direction]
  !! @param [out] Ny number  of components in the y direction
  !! @param [out] Ky number of singular vectors in the y direction
  !! @param [out] Nz number  of components in the z direction
  !! @param [out] Kz number of singular vectors in the z direction
  !!
  !! \details Negatives value means that there is no component in the given
  !!  direction.
  !<
  subroutine get_evd_cov_size(C, Nx, Kx, Ny, Ky, Nz, Kz)
    type(evd_cov_matrix), intent(in out) :: C
    INTEGER, INTENT(OUT), OPTIONAL :: Nx, Kx, Ny, Ky, Nz, Kz

    if(C%im_nFac>=1)then
      if( present(Nx) )Nx = size(C%Ux,1)
      if( present(Kx) )Kx = size(C%Ux,2)
    else
      if( present(Nx) )Nx = -1
      if( present(Kx) )Kx = -1
    end if

    if(C%im_nFac>=2)then
      if( present(Ny) )Ny = size(C%Uy,1)
      if( present(Ky) )Ky = size(C%Uy,2)
    else
      if( present(Ny) )Ny = -1
      if( present(Ky) )Ky = -1
    end if

    if(C%im_nFac>=3)then
      if( present(Nz) )Nz = size(C%Uz,1)
      if( present(Kz) )Kz = size(C%Uz,2)
    else
      if( present(Nz) )Nz = -1
      if( present(Kz) )Kz = -1
    end if
  end subroutine get_evd_cov_size

  !> \brief applies the covariance matrix to the vector f, g = Bf
  !! @param [in] C covariance matrix derived type
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behavior of the program if the covariance
  !!  matrix is not loaded. if TRUE, the program stops, otherwise g =f.
  !!  The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_evd_cov(C, f, g, die)
    type(evd_cov_matrix), intent(in) :: C
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables
    LOGICAL :: ll_die

    IF( PRESENT(die) )THEN
      ll_die = die
    ELSE
      ll_die = .TRUE.
    END IF

    SELECT CASE(C%im_nFac)
      CASE (1)
        CALL cov_mv(C%Ux, C%Dx, f, g)
      CASE (2)
        CALL cov_mv(C%Ux, C%Dx, C%Uy, C%Dy, f, g)
      CASE (3)
        CALL cov_mv(C%Ux, C%Dx, C%Uy, C%Dy, C%Uz, C%Dz, f, g)
      CASE (default_nFact)!default value
        IF(ll_die)THEN
          CALL stop_program('In apply_evd_cov covariance matrix is not loaded')
        ELSE
          g=f
        END IF
       case default
          call debug(C%im_nFac, "Number of factor: ")
          CALL stop_program('In apply_evd_cov covariance matrix is not supported')
    END SELECT
  END SUBROUTINE apply_evd_cov

  !> \brief applies the inverse of covariance matrix, g = B^-1 f
  !! @param [in] C covariance matrix derived type
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance
  !!   matrix is not loaded. if TRUE, the program stops, otherwise g =f.
  !!   The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_evd_cov_inv(C, f, g, die)
    type(evd_cov_matrix), intent(in) :: C
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables
    LOGICAL :: ll_die

    IF( PRESENT(die) )THEN
      ll_die = die
    ELSE
      ll_die = .TRUE.
    END IF

    SELECT CASE(C%im_nFac)
      CASE (1)
        CALL cov_inv_mv(C%Ux, C%Dx, f, g)
      CASE (2)
        CALL cov_inv_mv(C%Ux, C%Dx, C%Uy, C%Dy, f, g)
      CASE (3)
        CALL cov_inv_mv(C%Ux, C%Dx, C%Uy, C%Dy, C%Uz, C%Dz, f, g)
      CASE default!default value
        IF(ll_die)THEN
          CALL stop_program('In apply_evd_cov_inv covariance matrix is not loaded')
        ELSE
          g=f
        END IF
    END SELECT
  END SUBROUTINE apply_evd_cov_inv

  !> \brief apply the square root of the covariance matrix to the vector f and
  !!  result g
  !! @param [in] C covariance matrix derived type
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance
  !!   matrix is not loaded. if TRUE, the program stops, otherwise g =f,
  !!   the user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_evd_square_cov(C, f, g, die)
    type(evd_cov_matrix), intent(in) :: C
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables
    LOGICAL :: ll_die

    IF( PRESENT(die) )THEN
      ll_die = die
    ELSE
      ll_die = .TRUE.
    END IF

    SELECT CASE(C%im_nFac)
      CASE (1)
        CALL cov_square_mv(C%Ux, C%Dx, f, g)
      CASE (2)
        CALL cov_square_mv(C%Ux, C%Dx, C%Uy, C%Dy, f, g)
      CASE (3)
        CALL cov_square_mv(C%Ux, C%Dx, C%Uy, C%Dy, C%Uz, C%Dz, f, g)
      CASE default!default value
        IF(ll_die)THEN
          CALL stop_program('In apply_evd_square_cov, covariance matrix is not loaded/supported')
        ELSE
          g=f
        END IF
    END SELECT
  END SUBROUTINE apply_evd_square_cov

  !> \brief apply the adjoint of the square root of the covariance matrix to
  !!   the vector f and result g
  !! @param [in] C covariance matrix derived type
  !! @param [in] fb input vector
  !! @param [in] gb output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance
  !!   matrix is not loaded. if TRUE, the program stops, otherwise g =f.
  !!   The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_evd_square_covb(C, fb, gb, die)
    type(evd_cov_matrix), intent(in) :: C
    REAL(cp), DIMENSION(:), INTENT(INOUT)  :: fb
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: gb
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables
    LOGICAL :: ll_die

    IF( PRESENT(die) )THEN
      ll_die = die
    ELSE
      ll_die = .TRUE.
    END IF

    SELECT CASE(C%im_nFac)
      CASE (1)
        CALL cov_square_mvb(C%Ux, C%Dx, fb, gb)
      CASE (2)
        CALL cov_square_mvb(C%Ux, C%Dx, C%Uy, C%Dy, fb, gb)
      CASE default !default value
        IF(ll_die)THEN
          CALL stop_program('In apply_evd_square_covb, covariance matrix is not loaded/supported')
        ELSE
          fb=fb+gb
          gb=0.0_cp
        END IF
    END SELECT
  END SUBROUTINE apply_evd_square_covb

  !> \brief apply the inverse of the square root of the covariance matrix
  !!   to the vector f and result g
  !! @param [in] C covariance matrix derived type
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance
  !!   matrix is not loaded. if TRUE, the program stops, otherwise g =f.
  !!   The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_evd_square_cov_inv(C, f, g, die)
    type(evd_cov_matrix), intent(in) :: C
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables
    LOGICAL :: ll_die

    IF( PRESENT(die) )THEN
      ll_die = die
    ELSE
      ll_die = .TRUE.
    END IF

    SELECT CASE(C%im_nFac)
      CASE (1)
        CALL cov_square_inv_mv( C%Ux, C%Dx, f, g )
      CASE (2)
        CALL cov_square_inv_mv( C%Ux, C%Dx, C%Uy, C%Dy, f, g )
      CASE (3)
        CALL cov_square_inv_mv( C%Ux, C%Dx, C%Uy, C%Dy, C%Uz, C%Dz, f, g )
      CASE ( default_nFact )!default value
        IF (ll_die) THEN
          CALL stop_program('In apply_evd_square_cov_inv, covariance matrix is not loaded')
        ELSE
          g=f
        END IF
    END SELECT
  END SUBROUTINE apply_evd_square_cov_inv

  !> @brief Makes a EVD covariance matrix
  !! @param [in, out] C is an evd covariance matrix data structure
  !! @param [in] xcp covariance parameters in the x dimension
  !! @param [in] ycp covariance parameters in the y dimension
  !! @param [in] zcp covariance parameters in the z dimension
  !!
  !! @details zcp is processed only if ycp is present
  !<
  subroutine make_evd_cov(C, xcp, ycp, zcp )
    type(evd_cov_matrix), intent(in out) :: C
    type(exp_cov_param), intent(in) :: xcp
    type(exp_cov_param), optional, intent(in) :: ycp, zcp
    !local variables
    real(dp), dimension(:,:), allocatable :: A

    allocate( A(xcp%n, xcp%n), C%Ux(xcp%n, xcp%k), C%Dx(xcp%k) )
    call se_cov_matrix(A, xcp)
    !call debug(A, 'in  make_evd_cov, Cx = ')
    call symmetric_ev(A, xcp%k, C%Ux, C%Dx)
    deallocate(A)
    C%im_nFac = 1

    if(present(ycp))then
      allocate( A(ycp%n, ycp%n), C%Uy(ycp%n, ycp%k), C%Dy(ycp%k) )
      call se_cov_matrix(A, ycp)
      call symmetric_ev(A, ycp%k, C%Uy, C%Dy)
      deallocate(A)
      C%im_nFac = C%im_nFac + 1
      if(present(zcp))then
        allocate( A(zcp%n, zcp%n), C%Uz(zcp%n, zcp%k), C%Dz(zcp%k) )
        call se_cov_matrix(A, zcp)
        call symmetric_ev(A, zcp%k, C%Uz, C%Dz)
        deallocate(A)
        C%im_nFac = C%im_nFac + 1
      end if
    else
      if(present(zcp))then
        call stop_program("In make_evd_cov: z is present and y is absent")
      end if
    end if

  end subroutine make_evd_cov

!end of public procedures

!privates procedures
  !> \brief Checks the parameters for the application of the covariance matrix
  !!  or its inverse
  !! @see cov_mv_2Fac for the description of the parameters
  !<
  SUBROUTINE cov_mv_check_2Fac(Ux, Dx, Uy, Dy, f, g)
    real(cp), dimension(:,:), intent(in) :: Ux, Uy
    real(cp), dimension(:), intent(in)  :: Dx, Dy, f
    real(cp), dimension(:), intent(in) :: g
    !local variables
    INTEGER :: nx, ny, kx, ky, nf_predict

    nx = SIZE(Ux,1)
    ny = SIZE(Uy,1)
    kx = SIZE(Ux,2)
    ky = SIZE(Uy,2)
    nf_predict = nx*ny
    if( (SIZE(Dx)/=kx).OR.(SIZE(Dy)/=ky).OR.(SIZE(f)/=nf_predict).OR.(SIZE(g)/=nf_predict) )THEN
      call debug('Expecting data with following specification, based on the shape of Ux and Uy', tag=dALLWAYS)
      call debug(kx, 'Dx of size :', tag=dALLWAYS)
      call debug(ky, 'Dy of size :', tag=dALLWAYS)
      call debug(nf_predict, 'f and g of size :', tag=dALLWAYS)
      call debug('Get the following, check and try aigain', tag=dALLWAYS)
      call debug(SIZE(Dx), 'SIZE(Dx)', tag=dALLWAYS)
      call debug(SIZE(Dy), 'SIZE(Dy)', tag=dALLWAYS)
      call debug(SIZE(f), 'SIZE(f)', tag=dALLWAYS)
      call debug(SIZE(g), 'SIZE(g)', tag=dALLWAYS)
      call debug('For debuging purpose, here are the shapes of Ux and Uy', tag=dALLWAYS)
      call debug(SHAPE(Ux), 'SHAPE(Ux)', tag=dALLWAYS)
      call debug(SHAPE(Uy), 'SHAPE(Uy)', tag=dALLWAYS)
      call stop_program('Non consistant input for the application of the covariance matrix or its inverse')
    end if
  end subroutine cov_mv_check_2Fac
  !> \brief Checks the parameters for the application of the covariance matrix
  !!  or its inverse
  !! @see cov_mv_3Fac for the description of the parameters
  !<
  SUBROUTINE cov_mv_check_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g)
    real(cp), dimension(:,:), intent(in) :: Ux, Uy, Uz
    real(cp), dimension(:), intent(in)  :: Dx, Dy, Dz, f
    real(cp), dimension(:), intent(in) :: g
    !local variables
    INTEGER :: nx, ny, nz, kx, ky, kz, nf_predict

    nx = size(Ux,1)
    ny = size(Uy,1)
    nz = size(Uz,1)
    kx = size(Ux,2)
    ky = size(Uy,2)
    kz = size(Uz,2)
    nf_predict = nx*ny*nz
    if( (SIZE(Dx)/=kx).OR.(SIZE(Dy)/=ky).OR.(SIZE(Dz)/=kz)&
        .OR.(SIZE(f)/=nf_predict).OR.(SIZE(g)/=nf_predict) )THEN
      call debug('Expecting data with following specification, based on the shape of Ux and Uy', tag=dALLWAYS)
      call debug(kx, 'Dx of size :', tag=dALLWAYS)
      call debug(ky, 'Dy of size :', tag=dALLWAYS)
      call debug(kz, 'Dz of size :', tag=dALLWAYS)
      call debug(nf_predict, 'f and g of size :', tag=dALLWAYS)
      call debug('Get the following, check and try aigain', tag=dALLWAYS)
      call debug(SIZE(Dx), 'SIZE(Dx)', tag=dALLWAYS)
      call debug(SIZE(Dy), 'SIZE(Dy)', tag=dALLWAYS)
      call debug(SIZE(Dz), 'SIZE(Dz)', tag=dALLWAYS)
      call debug(SIZE(f), 'SIZE(f)', tag=dALLWAYS)
      call debug(SIZE(g), 'SIZE(g)', tag=dALLWAYS)
      call debug('For debuging purpose, here are the shapes of Ux and Uy', tag=dALLWAYS)
      call debug(SHAPE(Ux), 'SHAPE(Ux)', tag=dALLWAYS)
      call debug(SHAPE(Uy), 'SHAPE(Uy)', tag=dALLWAYS)
      call debug(SHAPE(Uz), 'SHAPE(Uz)', tag=dALLWAYS)
      call stop_program('Non consistant input for the application of the covariance matrix or its inverse')
    end if
  end subroutine cov_mv_check_3Fac

  !> @brief multiply a vector by the covariance matrix given as a Kronecker
  !!   product of two covariance Bx and By. Each one being given by its EVD.
  !!   Best suited for two physical dimensions problems where the covariance
  !!   factorizes.
  !! @param[in] Ux matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix Bx in the direction x, size (nx, kx)
  !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
  !! @param[in] Uy matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix By in the direction y, size (ny, ky)
  !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
  !! @param[in] f input vector of size nx*ny (size of the discretization grid)
  !! @param[out] g output vector of size nx*ny
  !! @details
  !! computes g = Bf
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the y-leading and the column major vectorization
  !<
  subroutine cov_mv_2Fac(Ux, Dx, Uy, Dy, f, g)
    real(cp), dimension(:,:), intent(in) :: Ux, Uy
    real(cp), dimension(:), intent(in)  :: Dx, Dy, f
    real(cp), dimension(:), intent(out) :: g
    !local variables
    real(cp), dimension( SIZE(Dx)*SIZE(Dy) ) :: t1, t2

    !checking the input parameters
    call cov_mv_check(Ux, Dx, Uy, Dy, f, g)

    call kron_matvec(Ux, Uy, f , t1, trans=.TRUE.)
    call kron_matvec(Dx, Dy, t1, t2)
    call kron_matvec(Ux, Uy, t2, g , trans=.FALSE.)
  end subroutine cov_mv_2Fac

  !> @brief multiply a vector by the covariance matrix given as a Kronecker
  !!   product of 3 covariance Bx, By, Bz. Each one being given by its EVD.
  !!   Best suited for two physical dimensions problems where the covariance
  !!   factorizes.
  !! @param[in] Ux matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix Bx in the direction x, size (nx, kx)
  !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
  !! @param[in] Uy matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix By in the direction y, size (ny, ky)
  !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
  !! @param[in] Uz matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix By in the direction y, size (nz, kz)
  !! @param[in] Dz vector of leading eigenvalues in the direction z, size kz
  !! @param[in] f input vector of size nx*ny*nz (size of the discretization grid)
  !! @param[out] g output vector of size nx*ny*nz
  !! @details
  !! computes g = Bf
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the y-leading and the column major vectorization
  !<
  subroutine cov_mv_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g)
    real(cp), dimension(:,:), intent(in) :: Ux, Uy, Uz
    real(cp), dimension(:), intent(in)  :: Dx, Dy, Dz, f
    real(cp), dimension(:), intent(out) :: g
    !local variables
    real(cp), dimension( SIZE(Dx)*SIZE(Dy)*SIZE(Dz) ) :: t1, t2

    !checking the input parameters
    call cov_mv_check(Ux, Dx, Uy, Dy, Uz, Dz, f, g)

    call kron_matvec(Ux, Uy, Uz, f , t1, trans=.TRUE.)
    call kron_matvec(Dx, Dy, Dz, t1, t2)
    call kron_matvec(Ux, Uy, Uz, t2, g , trans=.FALSE.)
  end subroutine cov_mv_3Fac

  !> \brief Checks the parameters for the application of the covariance matrix
  !!   or its inverse
  !! @see cov_mv_2Fac for the description of the parameters
  !<
  subroutine cov_mv_check_default(U, D, f, g)
    real(cp), dimension(:,:), intent(in) :: U
    real(cp), dimension(:), intent(in)  :: D, f
    real(cp), dimension(:), intent(out) :: g
    !local variables
    integer :: n, k

    n = SIZE(U,1)
    k = SIZE(U,2)
    if( (SIZE(D)/=k).OR.(SIZE(f)/=n).OR.(SIZE(g)/=n) )THEN
      call debug('cov_mv_check_default: Expecting data with following specification, based on the shape of U', tag=dALLWAYS)
      call debug(k, 'D of size :', tag=dALLWAYS)
      call debug(n, 'f and g of size :', tag=dALLWAYS)
      call debug('Get the following, check and try aigain', tag=dALLWAYS)
      call debug(SIZE(D), 'SIZE(D)', tag=dALLWAYS)
      call debug(SIZE(f), 'SIZE(f)', tag=dALLWAYS)
      call debug(SIZE(g), 'SIZE(g)', tag=dALLWAYS)
      call debug('For debuging purpose, here are the shapes of Ux and Uy', tag=dALLWAYS)
      call debug(SHAPE(U), 'SHAPE(U)', tag=dALLWAYS)
      call stop_program('Non consistant input for the application of the covariance matrix or its inverse')
    end if
  end subroutine cov_mv_check_default

  !> @brief multiply a vector by the covariance matrix given as its EVD
  !! @param[in] U matrix whose columns are the leading eigenvectors of the
  !!   covariance matrix B, size (n, k)
  !! @param[in] D vector of leading eigenvalues, size k
  !! @param[in] f input vector of size n (size of the discretization grid)
  !! @param[out] g output vector of size n
  !! @details
  !! computes g = Bf
  !<
  subroutine cov_mv_default(U, D, f, g)
    real(cp), dimension(:,:), intent(in) :: U
    real(cp), dimension(:), intent(in)  :: D, f
    real(cp), dimension(:), intent(out) :: g
    !local variables
    REAL(cp), DIMENSION( SIZE(D) ) :: t1, t2

    !checking the input parameters
    CALL cov_mv_check( U, D, f, g )

    t1 = MATMUL(f, U) !transpose
    t2 = D*t1
    g  = MATMUL(U, t2)
  end subroutine cov_mv_default

  !> \brief multiplies a vector by the inverse of the covariance matrix,
  !!   two physical dimensions problem
  !! @param[in] Ux matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix Bx in the direction x, size (nx, kx)
  !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
  !! @param[in] Uy matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix By in the direction y, size (ny, ky)
  !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
  !! @param[in] f input vector of size nx*ny (size of the discretization grid)
  !! @param[out] g output vector of size nx*ny
  !! \details
  !! computes g = B^-1 f
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the y-leading and the column major vectorization
  !<
  SUBROUTINE cov_inv_mv_2Fac(Ux, Dx, Uy, Dy, f, g)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: Ux, Uy
    REAL(cp), DIMENSION(:), INTENT(IN)  :: Dx, Dy, f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    !local variables
    REAL(cp), DIMENSION(SIZE(Dx)) :: DxInv
    REAL(cp), DIMENSION(SIZE(Dy)) :: DyInv
    REAL(cp), DIMENSION( SIZE(Dx)*SIZE(Dy) ) :: t1, t2

    !checking the input parameters
    CALL cov_mv_check(Ux, Dx, Uy, Dy, f, g)

    DxInv = 1.0/Dx
    DyInv = 1.0/Dy
    CALL kron_matvec(Ux, Uy, f , t1, trans=.TRUE.)
    CALL kron_matvec(DxInv, DyInv, t1, t2)
    CALL kron_matvec(Ux, Uy, t2, g , trans=.FALSE.)
  END SUBROUTINE cov_inv_mv_2Fac

  !> @brief multiplies a vector by the inverse of the covariance matrix,
  !!   3 physical dimensions problem
  !! @param[in] Ux matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix Bx in the direction x, size (nx, kx)
  !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
  !! @param[in] Uy matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix By in the direction y, size (ny, ky)
  !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
  !! @param[in] Uz matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix By in the direction y, size (nz, kz)
  !! @param[in] Dz vector of leading eigenvalues in the direction z, size kz
  !! @param[in] f input vector of size nx*ny*nz (size of the discretization grid)
  !! @param[out] g output vector of size nx*ny*nz
  !! @details
  !! computes g = Bf
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the y-leading and the column major vectorization
  !<
  SUBROUTINE cov_inv_mv_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: Ux, Uy, Uz
    REAL(cp), DIMENSION(:), INTENT(IN)  :: Dx, Dy, Dz, f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    !local variables
    REAL(cp), DIMENSION(SIZE(Dx)) :: DxInv
    REAL(cp), DIMENSION(SIZE(Dy)) :: DyInv
    REAL(cp), DIMENSION(SIZE(Dz)) :: DzInv
    REAL(cp), DIMENSION( SIZE(Dx)*SIZE(Dy)*SIZE(Dz) ) :: t1, t2

    !checking the input parameters
    CALL cov_mv_check(Ux, Dx, Uy, Dy, Uz, Dz, f, g)

    DxInv = 1.0/Dx
    DyInv = 1.0/Dy
    DzInv = 1.0/Dz
    CALL kron_matvec(Ux, Uy, Uz, f , t1, trans=.TRUE.)
    CALL kron_matvec(DxInv, DyInv, DzInv, t1, t2)
    CALL kron_matvec(Ux, Uy, Uz, t2, g , trans=.FALSE.)
  END SUBROUTINE cov_inv_mv_3Fac

  !> \brief multiplies a vector by the inverse of the covariance matrix given
  !!   as its EVD
  !! @param[in] U matrix whose columns are the leading eigenvalues vector of
  !!   the covariance matrix B, size (n, k)
  !! @param[in] D vector of leading eigenvalues, size k
  !! @param[in] f input vector of size n (size of the discretization grid)
  !! @param[out] g output vector of size n
  !! \details
  !! computes g = B^-1 f
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the y-leading and the column major vectorization
  !<
  SUBROUTINE cov_inv_mv_default(U, D, f, g)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: U
    REAL(cp), DIMENSION(:), INTENT(IN)  :: D, f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    !local variables
    REAL(cp), DIMENSION(SIZE(D)) :: t1, t2

    !checking the input parameters
    CALL cov_mv_check( U, D, f, g )

    t1 = MATMUL(f, U) !transpose
    t2 = t1/D
    g  = MATMUL(U, t2)
  END SUBROUTINE cov_inv_mv_default

  !> \brief Checks the parameters for the application of the square root of
  !!  the covariance matrix or its inverse
  !! @see cov_square_mv_default for the description of the parameters
  !<
  SUBROUTINE cov_square_mv_check_default(U, D, f, g, inv)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: U
    REAL(cp), DIMENSION(:), INTENT(IN)  :: D, f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL :: inv
    !local variables
    INTEGER :: n, k, nf_predict, ng_predict
    LOGICAL :: ll_inv

    IF( PRESENT(inv) )THEN
      ll_inv = inv
    ELSE
      ll_inv = .FALSE.
    END IF
    n = SIZE(U,1)
    k = SIZE(U,2)
    IF(ll_inv)THEN
      nf_predict = n
      ng_predict = k
    ELSE
      nf_predict = k
      ng_predict = n
    END IF
    IF( (SIZE(D)/=k).OR.(SIZE(f)/=nf_predict).OR.(SIZE(g)/=ng_predict) )THEN
      CALL debug('Expecting data with following specification, based on the shape of Ux and Uy', tag=dALLWAYS)
      CALL debug(k, 'D of size :', tag=dALLWAYS)
      CALL debug(nf_predict, 'f of size :', tag=dALLWAYS)
      CALL debug(ng_predict, 'g of size :', tag=dALLWAYS)
      CALL debug('Get the following, check and try aigain', tag=dALLWAYS)
      CALL debug(SIZE(D), 'SIZE(Dx)', tag=dALLWAYS)
      CALL debug(SIZE(f), 'SIZE(f)', tag=dALLWAYS)
      CALL debug(SIZE(g), 'SIZE(g)', tag=dALLWAYS)
      CALL debug('For debuging purpose, here are the shapes of Ux and Uy', tag=dALLWAYS)
      CALL debug(SHAPE(U), 'SHAPE(U)', tag=dALLWAYS)
      IF(ll_inv) CALL debug('The purpose is to compute apply the inverse square root', tag=dALLWAYS)
      CALL stop_program('Non consistant input for the application of the covariance matrix or its inverse')
    END IF
  END SUBROUTINE cov_square_mv_check_default

  !> \brief Checks the parameters for the application of the square root
  !!   of the covariance matrix or its inverse
  !! @see cov_square_mv_2Fac for the description of the parameters
  !<
  SUBROUTINE cov_square_mv_check_2Fac(Ux, Dx, Uy, Dy, f, g, inv)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: Ux, Uy
    REAL(cp), DIMENSION(:), INTENT(IN)  :: Dx, Dy, f
    REAL(cp), DIMENSION(:), INTENT(IN) :: g
    LOGICAL, OPTIONAL :: inv
    !local variables
    INTEGER :: nx, ny, kx, ky, nf_predict, ng_predict
    LOGICAL :: ll_inv

    IF( PRESENT(inv) )THEN
      ll_inv = inv
    ELSE
      ll_inv = .FALSE.
    END IF
    nx = SIZE(Ux,1)
    ny = SIZE(Uy,1)
    kx = SIZE(Ux,2)
    ky = SIZE(Uy,2)
    IF(ll_inv)THEN
      nf_predict = nx*ny
      ng_predict = kx*ky
    ELSE
      nf_predict = kx*ky
      ng_predict = nx*ny
    END IF
    IF( (SIZE(Dx)/=kx).OR.(SIZE(Dy)/=ky).OR.(SIZE(f)/=nf_predict).OR.(SIZE(g)/=ng_predict) )THEN
      CALL debug('Expecting data with following specification, based on the shape of Ux and Uy', tag=dALLWAYS)
      CALL debug(kx, 'Dx of size :', tag=dALLWAYS)
      CALL debug(ky, 'Dy of size :', tag=dALLWAYS)
      CALL debug(nf_predict, 'f of size :', tag=dALLWAYS)
      CALL debug(ng_predict, 'g of size :', tag=dALLWAYS)
      CALL debug('Get the following, check and try aigain', tag=dALLWAYS)
      CALL debug(SIZE(Dx), 'SIZE(Dx)', tag=dALLWAYS)
      CALL debug(SIZE(Dy), 'SIZE(Dy)', tag=dALLWAYS)
      CALL debug(SIZE(f), 'SIZE(f)', tag=dALLWAYS)
      CALL debug(SIZE(g), 'SIZE(g)', tag=dALLWAYS)
      CALL debug('For debuging purpose, here are the shapes of Ux and Uy', tag=dALLWAYS)
      CALL debug(SHAPE(Ux), 'SHAPE(Ux)', tag=dALLWAYS)
      CALL debug(SHAPE(Uy), 'SHAPE(Uy)', tag=dALLWAYS)
      IF(ll_inv) CALL debug('The purpose is to apply the inverse square root', tag=dALLWAYS)
      CALL stop_program('Non consistant input for the application of the covariance matrix or its inverse')
    END IF
  END SUBROUTINE cov_square_mv_check_2Fac


  !> \brief Checks the parameters for the application of the square root
  !!   of the covariance matrix or its inverse
  !! @see cov_square_mv_3Fac for the description of the parameters
  !<
  SUBROUTINE cov_square_mv_check_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g, inv)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: Ux, Uy, Uz
    REAL(cp), DIMENSION(:), INTENT(IN)  :: Dx, Dy, Dz, f
    REAL(cp), DIMENSION(:), INTENT(IN) :: g
    LOGICAL, OPTIONAL :: inv
    !local variables
    INTEGER :: nx, ny, nz, kx, ky, kz, nf_predict, ng_predict
    LOGICAL :: ll_inv

    IF( PRESENT(inv) )THEN
      ll_inv = inv
    ELSE
      ll_inv = .FALSE.
    END IF
    nx = SIZE(Ux,1)
    ny = SIZE(Uy,1)
    nz = SIZE(Uz,1)
    kx = SIZE(Ux,2)
    ky = SIZE(Uy,2)
    kz = SIZE(Uz,2)
    IF(ll_inv)THEN
      nf_predict = nx*ny*nz
      ng_predict = kx*ky*kz
    ELSE
      nf_predict = kx*ky*kz
      ng_predict = nx*ny*nz
    END IF
    IF( (SIZE(Dx)/=kx).OR.(SIZE(Dy)/=ky).OR.(SIZE(Dz)/=kz)&
        .OR.(SIZE(f)/=nf_predict).OR.(SIZE(g)/=ng_predict) )THEN
      CALL debug('Expecting data with following specification, based on the shape of Ux and Uy', tag=dALLWAYS)
      CALL debug(kx, 'Dx of size :', tag=dALLWAYS)
      CALL debug(ky, 'Dy of size :', tag=dALLWAYS)
      CALL debug(kz, 'Dz of size :', tag=dALLWAYS)
      CALL debug(nf_predict, 'f of size :', tag=dALLWAYS)
      CALL debug(ng_predict, 'g of size :', tag=dALLWAYS)
      CALL debug('Got the following, check and try aigain', tag=dALLWAYS)
      CALL debug(SIZE(Dx), 'SIZE(Dx)', tag=dALLWAYS)
      CALL debug(SIZE(Dy), 'SIZE(Dy)', tag=dALLWAYS)
      CALL debug(SIZE(Dz), 'SIZE(Dz)', tag=dALLWAYS)
      CALL debug(SIZE(f), 'SIZE(f)', tag=dALLWAYS)
      CALL debug(SIZE(g), 'SIZE(g)', tag=dALLWAYS)
      CALL debug('For debuging purpose, here are the shapes of Ux, Uy, Uz', tag=dALLWAYS)
      CALL debug(SHAPE(Ux), 'SHAPE(Ux)', tag=dALLWAYS)
      CALL debug(SHAPE(Uy), 'SHAPE(Uy)', tag=dALLWAYS)
      CALL debug(SHAPE(Uz), 'SHAPE(Uz)', tag=dALLWAYS)
      IF(ll_inv) CALL debug('The purpose is to apply the inverse square root', tag=dALLWAYS)
      CALL stop_program('Non consistant input for the application of the covariance matrix or its inverse')
    END IF
  END SUBROUTINE cov_square_mv_check_3Fac

    !> \brief multiplies a vector by the square root of the covariance matrix
    !!  given as its EVD
    !! @param[in] U matrix whose columns are the leading eigenvalues vector
    !!   of the covariance matrix B, size (n, k)
    !! @param[in] D vector of leading eigenvalues, size k
    !! @param[in] f input vector of size k (size of the discretization grid)
    !! @param[out] g output vector of size n
    !! \details
    !! computes g = B^0.5 f
    !<
    SUBROUTINE cov_square_mv_default(U, D, f, g)
        REAL(cp), DIMENSION(:,:), INTENT(IN) :: U
        REAL(cp), DIMENSION(:), INTENT(IN)  :: D, f
        REAL(cp), DIMENSION(:), INTENT(OUT) :: g
        !local variables
        REAL(cp), DIMENSION( SIZE(D) ) :: t

        !checking the input parameters
        CALL cov_square_mv_check(U, D, f, g, inv=.FALSE.)

        t = SQRT(D)*f
        g = MATMUL(U, t)
    END SUBROUTINE cov_square_mv_default

    !> \brief multiplies a vector by the square root of the covariance matrix,
    !!  two physical dimensions problem
    !! @param[in] Ux matrix whose columns are the leading eigenvalues vector
    !!   of the covariance matrix Bx in the direction x, size (nx, kx)
    !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
    !! @param[in] Uy matrix whose columns are the leading eigenvalues vector
    !!   of the covariance matrix By in the direction y, size (ny, ky)
    !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
    !! @param[in] f input vector of size kx*ky (size of the discretization grid)
    !! @param[out] g output vector of size nx*ny
    !! \details
    !! computes g = B^0.5 f
    !! this routine assumes either :
    !! - the x-leading and and row major vectorization or
    !! - the y-leading and the column major vectorization
    !<
    SUBROUTINE cov_square_mv_2Fac(Ux, Dx, Uy, Dy, f, g)
        REAL(cp), DIMENSION(:,:), INTENT(IN) :: Ux, Uy
        REAL(cp), DIMENSION(:), INTENT(IN)  :: Dx, Dy, f
        REAL(cp), DIMENSION(:), INTENT(OUT) :: g
        !local variables
        REAL(cp), DIMENSION(SIZE(Dx)) :: Dx05
        REAL(cp), DIMENSION(SIZE(Dy)) :: Dy05
        REAL(cp), DIMENSION( SIZE(Dx)*SIZE(Dy) ) :: t

        !checking the input parameters
        CALL cov_square_mv_check_2Fac(Ux, Dx, Uy, Dy, f, g, inv=.FALSE.)

        Dx05 = SQRT(Dx)
        Dy05 = SQRT(Dy)
    !     call debug (Dx(size(Dx)-4:),'Dx = ', tag=dALLWAYS)
    !     call debug (Dy(size(Dy)-4:),'Dy = ', tag=dALLWAYS)
    !     call debug (Dx05(size(Dx05)-4:),'Dx05 = ', tag=dALLWAYS)
    !     call debug (Dy05(size(Dy05)-4:),'Dy05 = ', tag=dALLWAYS)
    !
    !     call debug (Dx(1:5),'Dx = ', tag=dALLWAYS)
    !     call debug (Dy(1:5),'Dy = ', tag=dALLWAYS)
    !     call debug (Dx05(1:5),'Dx05 = ', tag=dALLWAYS)
    !     call debug (Dy05(1:5),'Dy05 = ', tag=dALLWAYS)

        CALL kron_matvec(Dx05, Dy05, f, t)
        !call debug (t(1:5),'t = ', tag=dALLWAYS)
        CALL kron_matvec(Ux, Uy, t, g , trans=.FALSE.)
        !call debug (g(1:5),'g = ', tag=dALLWAYS)
    END SUBROUTINE cov_square_mv_2Fac

  !> \brief multiplies a vector by the square root of the covariance matrix,
  !!  two physical dimensions problem
  !! @param[in] Ux matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix Bx in the direction x, size (nx, kx)
  !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
  !! @param[in] Uy matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix By in the direction y, size (ny, ky)
  !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
  !! @param[in] Uz matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix Bz in the direction z, size (nz, kz)
  !! @param[in] Dz vector of leading eigenvalues in the direction z, size kz
  !! @param[in] f input vector of size kx*ky*kz (size of the discretization grid)
  !! @param[out] g output vector of size nx*ny*nz
  !! \details
  !! computes g = B^0.5 f
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the z-leading and the column major vectorization
  !<
  subroutine cov_square_mv_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g)
    real(cp), dimension(:,:), intent(in) :: Ux, Uy, Uz
    real(cp), dimension(:), intent(in)  :: Dx, Dy, Dz, f
    real(cp), dimension(:), intent(out) :: g
    !local variables
    real(cp), dimension(size(dx))       :: Dx05
    real(cp), dimension(size(dy))       :: Dy05
    real(cp), dimension(size(dz))       :: Dz05
    real(cp), dimension( size(Dx)*size(Dy)*size(Dz) ) :: t

    !checking the input parameters
    call cov_square_mv_check_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g, inv=.false.)

    Dx05 = sqrt(dx)
    dy05 = sqrt(Dy)
    dz05 = sqrt(Dz)

    call kron_matvec(Dx05, Dy05, Dz05, f, t)
    !call debug (t(1:5),'t = ', tag=dALLWAYS)
    call kron_matvec(Ux, Uy, Uz, t, g , trans=.false.)
    !call debug (g(1:5),'g = ', tag=dALLWAYS)
  end subroutine cov_square_mv_3Fac

  !> adjoint subroutine, used mathmatical formulation of the adjoint with
  !! linear operators
  !<
  SUBROUTINE cov_square_mvb_default(U, D, fb, gb)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: U
    REAL(cp), DIMENSION(:), INTENT(IN)  :: D
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: gb, fb
    !local variables
    REAL(cp), DIMENSION( SIZE(D) ) :: tb

    !checking the input parameters
    CALL cov_square_mv_check(U, D, fb, gb, inv=.FALSE.)

    tb = tb + MATMUL(gb, U)
    gb = 0.0_cp
    fb = fb + SQRT(D)*tb
    tb = 0.0_cp
  END SUBROUTINE cov_square_mvb_default

  !> adjoint subroutine, used mathmatical formulation of the adjoint with
  !! linear operators
  !<
  SUBROUTINE cov_square_mvb_2Fac(Ux, Dx, Uy, Dy, fb, gb)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: Ux, Uy
    REAL(cp), DIMENSION(:), INTENT(IN)  :: Dx, Dy
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: gb, fb
    !local variables
    REAL(cp), DIMENSION(SIZE(Dx)) :: Dx05
    REAL(cp), DIMENSION(SIZE(Dy)) :: Dy05
    REAL(cp), DIMENSION( SIZE(Dx)*SIZE(Dy) ) :: tb
    REAL(cp), DIMENSION( SIZE(fb) ) :: tmp

    !checking the input parameters
    CALL cov_square_mv_check_2Fac(Ux, Dx, Uy, Dy, fb, gb, inv=.FALSE.)
    !initialization

    Dx05 = SQRT(Dx)
    Dy05 = SQRT(Dy)

    CALL kron_matvec( Ux, Uy, gb, tb , trans=.TRUE. )!/!\optimised, knowing that tb is zero
    gb = 0.0_cp
    CALL kron_matvec( Dx05, Dy05, tb, tmp )! No need for transpose with diagonal square matrices
    fb = fb + tmp
    tb = 0.0_cp
  END SUBROUTINE cov_square_mvb_2Fac

  !> \brief multiplies a vector by the inverse of the square root of
  !!   the covariance matrix, two physical dimensions problem
  !! @param[in] Ux matrix whose columns are the leading eigenvalues
  !!   vector of the covariance matrix Bx in the direction x, size (nx, kx)
  !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
  !! @param[in] Uy matrix whose columns are the leading eigenvalues
  !!   vector of the covariance matrix By in the direction y, size (ny, ky)
  !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
  !! @param[in] f input vector of size nx*ny (size of the discretization grid)
  !! @param[out] g output vector of size nx*ny
  !! \details
  !! computes g = B^-0.5 f
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the y-leading and the column major vectorization
  !<
  SUBROUTINE cov_square_inv_mv_2Fac(Ux, Dx, Uy, Dy, f, g)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: Ux, Uy
    REAL(cp), DIMENSION(:), INTENT(IN)  :: Dx, Dy, f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    !local variables
    REAL(cp), DIMENSION(SIZE(Dx)) :: DxInv
    REAL(cp), DIMENSION(SIZE(Dy)) :: DyInv
    REAL(cp), DIMENSION( SIZE(Dx)*SIZE(Dy) ) :: t

    !checking the input parameters
    CALL cov_square_mv_check_2Fac(Ux, Dx, Uy, Dy, f, g, inv=.TRUE.)

    DxInv = SQRT(1.0/Dx)
    DyInv = SQRT(1.0/Dy)
    CALL kron_matvec(Ux, Uy, f, t , trans=.TRUE.)
    CALL kron_matvec(DxInv, DyInv, t, g)
  END SUBROUTINE cov_square_inv_mv_2Fac

  !> \brief multiplies a vector by the inverse of the square root of
  !!   the covariance matrix, 3 physical dimensions problem
  !! @param[in] Ux matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix Bx in the direction x, size (nx, kx)
  !! @param[in] Dx vector of leading eigenvalues in the direction x, size kx
  !! @param[in] Uy matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix By in the direction y, size (ny, ky)
  !! @param[in] Dy vector of leading eigenvalues in the direction y, size ky
  !! @param[in] Uz matrix whose columns are the leading eigenvalues vector
  !!   of the covariance matrix Bz in the direction z, size (nz, kz)
  !! @param[in] Dz vector of leading eigenvalues in the direction z, size kz
  !! @param[in] f input vector of size kx*ky*kz (size of the discretization grid)
  !! @param[out] g output vector of size nx*ny*nz
  !! \details
  !! computes g = B^-0.5 f
  !! this routine assumes either :
  !! - the x-leading and and row major vectorization or
  !! - the z-leading and the column major vectorization
  !<
  subroutine cov_square_inv_mv_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g)
    real(cp), dimension(:,:), intent(in) :: Ux, Uy, Uz
    real(cp), dimension(:), intent(in)  :: Dx, Dy, Dz, f
    real(cp), dimension(:), intent(out) :: g
    !local variables
    real(cp), dimension(size(dx))       :: DxInv
    real(cp), dimension(size(dy))       :: DyInv
    real(cp), dimension(size(dz))       :: DzInv
    real(cp), dimension( size(Dx)*size(Dy)*size(Dz) ) :: t

    !checking the input parameters
    call cov_square_mv_check_3Fac(Ux, Dx, Uy, Dy, Uz, Dz, f, g, inv=.false.)

    DxInv = SQRT(1.0/Dx)
    DyInv = SQRT(1.0/Dy)
    DzInv = SQRT(1.0/Dz)
    CALL kron_matvec(Ux, Uy, Uz, f, t , trans=.TRUE.)
    CALL kron_matvec(DxInv, DyInv, DzInv, t, g)
  end subroutine cov_square_inv_mv_3Fac

  !> \brief multiplies a vector by the inverse of the square root of
  !!   the covariance matrix given as its EVD
  !! @param[in] U matrix whose columns are the leading eigenvalues
  !!   vector of the covariance matrix B, size (n, k)
  !! @param[in] D vector of leading eigenvalues, size k
  !! @param[in] f input vector of size n (size of the discretization grid)
  !! @param[out] g output vector of size n
  !! \details
  !! computes g = B^-0.5 f
  !<
  SUBROUTINE cov_square_inv_mv_default(U, D, f, g)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: U
    REAL(cp), DIMENSION(:), INTENT(IN)  :: D, f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    !local variables
    REAL(cp), DIMENSION( SIZE(D) ) :: t

    !checking the input parameters
    CALL cov_square_mv_check(U, D, f, g, inv=.TRUE.)

    t = MATMUL(f, U)
    g = t/SQRT(D)
  END SUBROUTINE cov_square_inv_mv_default

  ! ***************************************************
  ! **************** debugging routines ***************
  ! ***************************************************

  !> \brief Check the conformity of the shape of matrices forming an EVD
  !! the program stops in case of non conformity
  !! @param [in] nU size of each singular vector
  !! @param [in] kU number singular vectors
  !! @param [in] kD1 number of singular values or on size of the diagonal matrix
  !! @param [in] kD2 optional (default 1) other size of the diagonal matrix
  !!   kD1 x kD2 is the shape of the diagonal matrix of singular values.
  !!   It can be stored as a full matrix or as a vector.
  !<
  SUBROUTINE EVD_check(nU, kU, kD1, kD2)
    INTEGER, INTENT(IN) :: nU, kU, kD1
    INTEGER, OPTIONAL, INTENT(IN) :: kD2
    !local variables
    INTEGER :: il_kD2

    IF( PRESENT(kD2)  )THEN
      il_kD2 = kD2
    ELSE
      il_kD2 = 1
    END IF

    IF( ( (kD1 == il_kD2).OR.( MIN(kD1, il_kD2) == 1 ) ) &
      .AND.( MAX(kD1, il_kD2) == kU ).AND.( nU >= kU ) )THEN
      !good configuration, nothing to do
    ELSE
      !bad configuration, message and stop the program
      CALL debug( 'Expecting EVD UDU^t to have the following specifications', tag=dALLWAYS )
      CALL debug( ' - U: number of rows greater or equal to the number of columns', tag=dALLWAYS )
      CALL debug( ' - D: vector or square matrix', tag=dALLWAYS )
      CALL debug( ' - MIN size of U equal to the MAX size of D', tag=dALLWAYS )
      CALL debug( 'For debuging purpose, here are the shapes of U and D', tag=dALLWAYS )
      CALL debug( (/nU, kU/), ' - SHAPE(U) = ', tag=dALLWAYS )
      CALL debug( (/kD1, il_kD2/), ' - SHAPE(D) = ', tag=dALLWAYS )
      CALL stop_program( 'Non consistant EVD' )
    END IF
  END SUBROUTINE EVD_check
!end of privates procedures
END MODULE evd_covariance
