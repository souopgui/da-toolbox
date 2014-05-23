!> \file kronecker.f90
!! |brief Module defining Kronecker operations for matrix computation
!!
!! @author Innocent Souopgui
!!
!<

!> Module defining Kronecker product for matrix computation
!!
!<
module kronecker
  USE general_constant
  USE debug_tools

  private kron_product_check, kron_product_check_diag, kron_product_check_full
  !>Forms the Kronecker product of two matrices
  INTERFACE kron_product
    MODULE PROCEDURE kron_product_full
    MODULE PROCEDURE kron_product_diag
  END INTERFACE kron_product

  !>Check the input for the Kronecker product of two matrix
  INTERFACE kron_product_check
    MODULE PROCEDURE kron_product_check_full
    MODULE PROCEDURE kron_product_check_diag
  END INTERFACE kron_product_check

  !> matrix vector multiplication with the matrix given by a Kronecker product of two matrices
  INTERFACE kron_matvec
    MODULE PROCEDURE kron_matvec_full ! matrices are full and represented by 2d arrays
    MODULE PROCEDURE kron_matvec_diag ! matrices are diagonal and represented by vectors
  END INTERFACE kron_matvec

CONTAINS

  !> \brief Computes the matrix vector product with the matrix given by a Kronecker product of two matrices.
  !! param[in] A first matrix of the kronecker product, size nrA x ncA
  !! param[in] B second matrix of the kronecker product, size nrB x ncB
  !! param[in] f input vector of size ncA*ncB
  !! param[out] g output vector of size nrA*nrB
  !! param[in] trans (optional) says if the transpose of the input matrices should be used instead of those matrices themselves.
  !! computes g = Kron(A,B)f if trans is not present or is set to TRUE, g = Kron(A,B)^Tf = (f^TKron(A,B) )^T if trans is present and set to TRUE
  !<
  SUBROUTINE kron_matvec_full(A, B, f, g, trans)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: A, B
    REAL(cp), DIMENSION(:), INTENT(IN)  :: f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    LOGICAL , INTENT(IN), OPTIONAL :: trans
    !local variables
    REAL(cp), DIMENSION(SIZE(B,1)) :: td!direct computation
    REAL(cp), DIMENSION(SIZE(B,2)) :: tt!transpose computation
    INTEGER :: nrA, ncA, nrB, ncB!number of row of X(nrX), !number of column of X(ncX)
    INTEGER :: iA, jA, if_start, if_end, ig_start, ig_end, nf_predicted, ng_predicted
    LOGICAL ll_trans

    IF( PRESENT(trans) )THEN
      ll_trans = trans
    ELSE
      ll_trans = .FALSE.
    END IF

    IF(ll_trans)THEN
      nf_predicted = SIZE(A,1)*SIZE(B,1)
      ng_predicted = SIZE(A,2)*SIZE(B,2)
    ELSE
      nf_predicted = SIZE(A,2)*SIZE(B,2)
      ng_predicted = SIZE(A,1)*SIZE(B,1)
    END IF
    IF( (SIZE(f)/=nf_predicted).OR.(SIZE(g)/=ng_predicted) )THEN
      CALL debug('Expecting data with following specification, based on the shape of A and B', tag=dALLWAYS)
      CALL debug(nf_predicted, 'f of size :', tag=dALLWAYS)
      CALL debug(ng_predicted, 'g of size :', tag=dALLWAYS)
      CALL debug('For debuging purpose, here are the shapes of Ux and Uy', tag=dALLWAYS)
      CALL debug(SHAPE(A), 'SHAPE(A)', tag=dALLWAYS)
      CALL debug(SHAPE(B), 'SHAPE(B)', tag=dALLWAYS)
      CALL stop_program('stopped from kron_matvec (with full matrices)')
    END IF

    nrA = SIZE(A,1)
    ncA = SIZE(A,2)
    nrB = SIZE(B,1)
    ncB = SIZE(B,2)

    g=0.0_cp
    IF(ll_trans)THEN! used the transposed
      if_start = 1
      if_end   = nrB
      DO iA = 1, nrA
        tt = MATMUL(f(if_start:if_end), B)
        ig_start = 1
        ig_end   = ncB
        DO jA=1,ncA
          g(ig_start:ig_end) = g(ig_start:ig_end) + A(iA,jA)*tt
          ig_start = ig_start + ncB
          ig_end   = ig_end   + ncB
        END DO
        if_start = if_start + nrB
        if_end   = if_end   + nrB
      END DO
    ELSE
      if_start = 1
      if_end   = ncB
      DO jA = 1, ncA
        td = MATMUL(B, f(if_start:if_end))
        ig_start = 1
        ig_end   = nrB
        DO iA=1,nrA
          g(ig_start:ig_end) = g(ig_start:ig_end) + A(iA,jA)*td
          ig_start = ig_start + nrB
          ig_end   = ig_end   + nrB
        END DO
        if_start = if_start + ncB
        if_end   = if_end   + ncB
      END DO
    END IF
  END SUBROUTINE kron_matvec_full

  !> \brief Computes the matrix vector product with the matrix given by a Kronecker product of two matrices
  !! the difference with kron_matvec_full is that the matrices are diagonal and represented by vectors of diagonal elements.
  !! @see kron_matvec_full for the description of parameters
  !<
  SUBROUTINE kron_matvec_diag(A, B, f, g)
    REAL(cp), DIMENSION(:), INTENT(IN) :: A, B, f
    REAL(cp), DIMENSION(:), INTENT(OUT) :: g
    !local variables
    INTEGER :: nA, nB !number of elements of X
    INTEGER :: iA, i_start, i_end, n_predicted


    nA = SIZE(A)
    nB = SIZE(B)
    n_predicted = nA*nB
    IF( (SIZE(f)/=n_predicted).OR.(SIZE(g)/=n_predicted) )THEN
      CALL debug('Expecting data with following specification, based on the shape of A and B', tag=dALLWAYS)
      CALL debug(n_predicted, 'f and g of size :', tag=dALLWAYS)
      CALL debug('For debuging purpose, here are the SIZE of A and B', tag=dALLWAYS)
      CALL debug(SIZE(A), 'SIZE(A)', tag=dALLWAYS)
      CALL debug(SIZE(B), 'SIZE(B)', tag=dALLWAYS)
      CALL stop_program('stopped from kron_matvec_diag')
    END IF

    g=0.0_cp

    i_start = 1
    i_end   = nB
    DO iA = 1, nA
      g(i_start:i_end) = A(iA)*( B*f(i_start:i_end) )
      i_start = i_start + nB
      i_end   = i_end   + nB
    END DO
  END SUBROUTINE kron_matvec_diag

  !> Checks parameters of the kronecker product
  !! @see kron_product_full
  SUBROUTINE kron_product_check_full(A, B, C)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: A, B, C
    !local variables
    INTEGER :: nrA, ncA, nrB, ncB!number of row of X(nrX), !number of column of X(ncX)

    nrA = SIZE(A,1)
    ncA = SIZE(A,2)
    nrB = SIZE(B,1)
    ncB = SIZE(B,2)
    IF( ( SIZE(C,1)/=nrA*nrB ).OR.( ( SIZE(C,2)/=ncA*ncB ) ) )THEN
      CALL debug('Expecting C with following specification, based on the shape of A and B', tag=dALLWAYS)
      CALL debug((/nrA*nrB, ncA*ncB/), 'expected SHAPE(C)', tag=dALLWAYS)
      CALL debug('Get the following, check and try aigain', tag=dALLWAYS)
      CALL debug(SHAPE(C), 'SHAPE(C)', tag=dALLWAYS)
      CALL debug('For debuging purpose, here are the shapes of A and B', tag=dALLWAYS)
      CALL debug(SHAPE(A), 'SHAPE(A)', tag=dALLWAYS)
      CALL debug(SHAPE(B), 'SHAPE(B)', tag=dALLWAYS)
      CALL stop_program('Bad shape for output matrix in kronecker product')
    END IF
  END SUBROUTINE kron_product_check_full

  !> Checks parameters of the kronecker product
  !! @see kron_product_full
  SUBROUTINE kron_product_check_diag(A, B, C)
    REAL(cp), DIMENSION(:), INTENT(IN) :: A, B, C
    !local variables
    INTEGER :: nA, nB

    nA = SIZE(A)
    nB = SIZE(B)
    IF( SIZE(C)/=nA*nB )THEN
      CALL debug('Expecting C with following specification, based on the SIZE of A and B', tag=dALLWAYS)
      CALL debug(nA*nB, 'expected SIZE(C)', tag=dALLWAYS)
      CALL debug('Get the following, check and try aigain', tag=dALLWAYS)
      CALL debug(SIZE(C), 'SIZE(C)', tag=dALLWAYS)
      CALL debug('For debuging purpose, here are the SIZEs of A and B', tag=dALLWAYS)
      CALL debug(SIZE(A), 'SIZE(A)', tag=dALLWAYS)
      CALL debug(SIZE(B), 'SIZE(B)', tag=dALLWAYS)
      CALL stop_program('Bad SIZE for output matrix in diagonal kronecker product')
    END IF
  END SUBROUTINE kron_product_check_diag

  !> \brief Forms the Kronecker product of two matrices
  !! \param[in] A first matrix in the Kronecker product
  !! \param[in] B second matrix in the Kronecker product
  !! \param[out] C resulting matrix
  !<
  SUBROUTINE kron_product_full(A, B, C)
    REAL(cp), DIMENSION(:,:), INTENT(IN) :: A, B
    REAL(cp), DIMENSION(:,:), INTENT(OUT) :: C
    !local variables
    INTEGER :: nrA, ncA, nrB, ncB, nrC, ncC!number of row of X(nrX), !number of column of X(ncX)
    INTEGER :: iA, jA, jB, iC_start, iC_end, jC

    !checking the input parameters
    CALL kron_product_check(A, B, C)
    nrA = SIZE(A,1)
    ncA = SIZE(A,2)
    nrB = SIZE(B,1)
    ncB = SIZE(B,2)
    nrC = SIZE(C,1)
    ncC = SIZE(C,2)

    DO jA=1,ncA
      !computing a bloc of columns of C
      jC = (jA-1)*ncB + 1
      DO jB=1,ncB
        !computing a column of C
        iC_start = 1
        iC_end   = nrB
        DO iA=1,nrA
          C(ic_start:ic_end, jC) = A(iA,jA)*B(:,jB)
          iC_start = iC_start + nrB
          iC_end   = iC_end   + nrB
        END DO
        jC = jC + 1
      END DO
    END DO
  END SUBROUTINE kron_product_full

  !> \brief Forms the Kronecker product of two matrices
  !! the difference with kron_matvec_full is that the matrices are diagonal and represented by vectors of diagonal elements.
  !! @see kron_product_full for the description of parameters
  !<
  SUBROUTINE kron_product_diag(A, B, C)
    REAL(cp), DIMENSION(:), INTENT(IN) :: A, B
    REAL(cp), DIMENSION(:), INTENT(OUT) :: C
    !local variables
    INTEGER :: nA, nB
    INTEGER :: iA, iC_start, iC_end

    !checking the input parameters
    CALL kron_product_check(A, B, C)
    nA = SIZE(A)
    nB = SIZE(B)

    iC_start = 1
    iC_end   = nB
    DO iA=1,nA
      C(ic_start:ic_end) = A(iA)*B
      iC_start = iC_start + nB
      iC_end   = iC_end   + nB
    END DO
  END SUBROUTINE kron_product_diag

end module kronecker