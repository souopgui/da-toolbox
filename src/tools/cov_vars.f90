!> \file cov_vars.f90
!! @brief  Module for covariance matrix variables; @see covariance.f90 for the decomposition of the covariance matrix
!! This module defines the variables representing the covariance matrix and the interface to apply the associated operators to a vector
!! It is assumed that the covariance matrix is provided as its Singular Value Decomposition: B = UDUt. Furthermore, it can be (Kronecker) factorized along the physical dimensions of the problem: B = Bx kron By kron Bz
!! @author Innocent Souopgui
!!
!! \todo rewrite programs/modules that use this module, to use directly the covariance module and define a local covariance matrix instead of the global one defined by this module
!!
MODULE cov_vars
  USE general_constant
  USE general_tools
  USE evd_covariance
IMPLICIT NONE
  PRIVATE
!   !> @brief  default number of factors in the Kronecker factorization of the covariance matrix, used for non initialized covariance matrix
!   !<
!   INTEGER, PARAMETER :: default_nFact=-1
!
!   !> @brief  number of terms in the Kronecker factorization of the covariance matrix.
!   !!Conceptually, it corresponds to the number of physical dimensions of the problem if the covariance is factorized.
!   !!Logically, it can be any Kronecker factorization.
!   !<
!   INTEGER, SAVE :: im_nFac = default_nFact
!
!   !> @brief  matrix(ces) of singular vector of the covariance matrix
!   !! If the covariance matrix B is not factorized, only Ux is initialized
!   !! IF B is factorized as B = Bx kron By (two terms), Ux and Uy are provided
!   !! A maximum of 3 terms in the factorization is supported. Only 1 and 2 terms are implemented for now.
!   !<
!   REAL(cp), ALLOCATABLE, DIMENSION(:,:), SAVE :: Ux, Uy, Uz
!
!   !> vector(s) of singular values of the EVD of B
!   REAL(cp), ALLOCATABLE, DIMENSION(:), SAVE :: Dx, Dy, Dz

  !> covariance matrix derived type
  type(evd_cov_matrix), save:: C

  PUBLIC load_covariance, apply_cov, apply_cov_inv, apply_square_cov, apply_square_covb, apply_square_cov_inv, cov_nsv

CONTAINS

  !> @brief  returns the number of singular values used to define the approximation of the covariance matrix; that is also the size of the changed state variable
  !! @param [in] [nctl] size of the control vector
  !! If a call to load_covariance has not been done before: the function returns nctl if it is provided, otherwise, the program is stopped. If a call to load_covariance have been done before, the associated value is returned
  !<
  FUNCTION cov_nsv(nctl) RESULT(nsv)
    INTEGER, OPTIONAL, INTENT(IN) :: nctl
    !local variables
    INTEGER :: nsv

    nsv = get_evd_cov_nsv(C, nctl)
  END FUNCTION cov_nsv

  !> @brief  Loads the covariance matrix represented by its singular value decomposition
  !! @param [in] dir directory containing the data files
  !! @param [out] Nx number of components [in the x direction]
  !! @param [out] Kx number of singular vectors [in the x direction]
  !! @param [out] Ny number  of components in the y direction
  !! @param [out] Ky number of singular vectors in the y direction
  !!
  !! It is assumed that the covariance matrix B is provided as its Singular Value Decomposition: B = UDU^t. Furthermore, it can be (Kronecker) factorized (for example along the physical dimensions of the problem) in which case it is given under the form:
  !! - B = Bx kron By = (UxDxUx^t) kron (UyDyUy^t) or
  !! - B = Bx kron By kron Bz = (UxDxUx^t) kron (UyDyUy^t) kron (UzDzUz^t)
  !! Conventions for the file names. Data files are named after the matrices names with the extension 'dat'
  !! - 'U.dat' and 'D.dat' for the non factorized case
  !! - 'Ux.dat',  'Uy.dat', 'Uz.dat' and 'Dx.dat', 'Dy.dat', 'Dz.dat' for the factorized case
  !! the existence of files is checked to determine the appropriate configuration.
  !! ('U.dat' and 'D.dat') or ('Ux.dat' and 'Dx.dat') should be present, otherwise the program terminates
  !!
  !<
  SUBROUTINE load_covariance(dir, Nx, Kx, Ny, Ky)
    CHARACTER(len=*), INTENT(IN)  :: dir
    INTEGER, INTENT(OUT), OPTIONAL :: Nx, Kx, Ny, Ky
    !local variables

    call load_evd_covariance(C, dir)
    call get_evd_cov_size(C, Nx, Kx, Ny, Ky)
  END SUBROUTINE load_covariance

  !> @brief  applies the covariance matrix to the vector f, g = Bf
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance matrix is not loaded. if TRUE, the program stops, otherwise g =f. The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_cov(f, g, die)
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables

    call apply_evd_cov(C, f, g, die)
  END SUBROUTINE apply_cov

  !> @brief  applies the inverse of covariance matrix, g = B^-1 f
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance matrix is not loaded. if TRUE, the program stops, otherwise g =f. The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_cov_inv(f, g, die)
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables

    call apply_evd_cov_inv(C, f, g, die)
  END SUBROUTINE apply_cov_inv

  !> @brief  apply the square root of the covariance matrix to the vector f and result g
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance matrix is not loaded. if TRUE, the program stops, otherwise g =f, the user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_square_cov(f, g, die)
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables

    call apply_evd_square_cov(C, f, g, die)
  END SUBROUTINE apply_square_cov

  !> @brief  apply the adjoint of the square root of the covariance matrix to the vector f and result g
  !! @param [in] fb input vector
  !! @param [in] gb output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance matrix is not loaded. if TRUE, the program stops, otherwise g =f. The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_square_covb(fb, gb, die)
    REAL(cp), DIMENSION(:), INTENT(INOUT)  :: fb
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: gb
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables

    call apply_evd_square_covb(C, fb, gb, die)
  END SUBROUTINE apply_square_covb

  !> @brief  apply the inverse of the square root of the covariance matrix to the vector f and result g
  !! @param [in] f input vector
  !! @param [in] g output vector
  !! @param [in] [die] specify the behaviour of the program if the covariance matrix is not loaded. if TRUE, the program stops, otherwise g =f. The user should make sure that f and g have the same size.
  !<
  SUBROUTINE apply_square_cov_inv(f, g, die)
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL, OPTIONAL, INTENT(IN) :: die
    !local variables

    call apply_evd_square_cov_inv(C, f, g, die)
  END SUBROUTINE apply_square_cov_inv

END MODULE cov_vars