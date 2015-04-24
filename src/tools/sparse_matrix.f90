!> @file sparse_matrix.f90
!! @brief Sparse matrix module
!! @author Innocent Souopgui
!! @details
!! The sparse matrix module provide the sparse matrix data structure
!! and some basic operation
!!
!! For a complete sparse matrix toolbox,
!! http://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html
!! Most subroutines in this module are borrowed from sparsekit
!<

!> @brief Module sparse matrix data structures and some basic subroutines
!!
!<
module sparse_matrix
    use general_constant
    use debug_tools
implicit none

    !> @brief user defined type for compressed sparse row matrix
    !! Let @a A be the matrix to store
    !! The @a val vector stores the values of the nonzero elements of
    !!   the matrix @a A as they are traversed in a row-wise fashion
    !! The @a col_ind vector stores the column indexes of the elements
    !!   in the @a val vector. If @a val[k] = @a A[i,j], then col_ind[k] = j
    !! The @a row_ptr vector stores the locations in the @a val vector that
    !!   start a row. If @a val[k] = @a A[i,j],
    !!   then row_ptr[i] <= k < row_ptr[i+1].
    !!   By convention, row_ptr[nrow+1] = nnz+1, where nnz is the number
    !!   of nonzeros in the matrix @a A
    !!
    !! @details http://web.eecs.utk.edu/~dongarra/etemplates/node373.html
    !<
    type CSR_matrix
        integer :: nrow = 0 !< Number of rows of the matrix
        integer :: ncol = 0 !< Number of columns of the matrix
        integer :: nnz = 0 !< Number of nonzero elements of the matrix
        integer, dimension(:), pointer :: col_ind => null()
        integer, dimension(:), pointer :: row_ptr => null()
        real(dp), dimension(:), pointer :: val => null()
        !> \brief indicates if the structure as been initialized
        logical :: initialized      = .false.
    end type CSR_matrix

contains

    !> @brief create a CSR matrix from a full matrix and a mask
    !! @param [in] A full matrix to compress
    !! @param [in] M mask matrix indicating masked entries of @a A
    !! @param [in, out] B csr matrix
    !! @details
    !!   @a A and @a M must have the same shape
    !<
    subroutine full2csr(A, M, B)
        real(dp), dimension(:,:), intent(in) :: A
        integer, dimension(:,:), intent(in)  :: M
        type(CSR_matrix), intent(inout) :: B
        !local variables
        integer :: nrow, ncol, nnz, i, j, pos

        if( sum( abs( shape(A)-shape(M) ) )/=0 )then
            call debug(shape(A), 'Shape of the input matrix')
            call debug(shape(M), 'Shape of the input mask  ')
            call stop_program('in to_csr, input matrix and mask have different shapes')
        end if
        nrow = size(A,1)
        ncol = size(A,2)
        nnz  = count( M == UNMASKED )
        call set_csr(B, nrow, ncol, nnz)
        if(nnz>0)then
            pos = 1!position in the val vector
            do i = 1, nrow
                B%row_ptr(i) = pos
                do j = 1, ncol
                    if(M(i,j)==UNMASKED)then
                        B%val(pos) = A(i,j)
                        B%col_ind(pos) = j
                        pos = pos+1
                    end if
                end do
            end do
            B%row_ptr(nrow+1) = nnz+1
        end if
    end subroutine full2csr

    !> @brief convert a CSR matrix to a full matrix
    !! @param [in] A csr matrix to expand
    !! @param [out] B full matrix (output)
    !! @param [out] M (optional) mask matrix indicating masked entries of @a B
    !! @details
    !!   @a A and @a M must have the same shape that corresponds
    !!   to the shape of the CSR matrix
    !<
    subroutine csr2full(A,B,M)
        type(CSR_matrix), intent(inout) :: A
        real(dp), dimension(:,:), intent(in out) :: B
        integer, dimension(:,:), optional, intent(in out)  :: M

        !local variables
        integer :: nrow, ncol, nnz, i, j, k

        nrow = A%nrow
        ncol = A%ncol
        if( sum( abs( shape(B)-[nrow,ncol] ) )/=0 )then
            call debug([A%nrow,A%ncol], 'Shape of the CSR matrix')
            call debug(shape(B), 'Shape of the full matrix  ')
            call stop_program('In to_csr, shapes of matrices differ')
        end if

        if( present(M) )then
            if( sum( abs( shape(M)-[nrow,ncol] ) )/=0 )then
                call debug([A%nrow,A%ncol], 'Shape of the CSR matrix')
                call debug(shape(M), 'Shape of the mask matrix  ')
                call stop_program('In to_csr, shapes of matrices differ')
            end if

            M = MASKED !suppose all masked
            do i = 1, A%nrow
                do k = A%row_ptr(i), A%row_ptr(i+1)-1
                    M(i, A%col_ind(k) ) = UNMASKED
                end do
            end do
        end if

        B = 0.0_dp
        do i = 1, A%nrow
            do k = A%row_ptr(i), A%row_ptr(i+1)-1
                B(i, A%col_ind(k) ) = A%val(k)
            end do
        end do
    end subroutine csr2full


    !> @brief create a CSR matrix from a full matrix and a mask
    !! @param [in, out] A sparse matrix data structure
    !! @param [in] nrow number of row of the matrix
    !! @param [in] ncol number of column of the matrix
    !! @param [in] nnz number of nonzero entry of the matrix
    !<
    subroutine new_csr(A, nrow, ncol, nnz)
        type(CSR_matrix), intent(inout) :: A
        integer, intent(in)             :: nrow, ncol, nnz

        call set_csr(A, nrow, ncol, nnz)
    end subroutine new_csr

    !> @brief Allocates state dynamic space for the csr matrix
    !! @param [in, out] A sparse matrix data structure
    !! @param [in] nrow number of row of the matrix
    !! @param [in] ncol number of column of the matrix
    !! @param [in] nnz number of nonzero entry of the matrix
    !! @details
    !! This routine is used to allocate space compressed sparse matrices.
    !<
    subroutine set_csr(A, nrow, ncol, nnz)
        type(CSR_matrix), intent(inout) :: A
        integer, intent(in)             :: nrow, ncol, nnz
        !local variables
        logical :: ll_allocate, ll_reset

        ll_allocate = .false.
        ll_reset = .false.

        if( (nnz/=0).and.( (nrow==0).or.(ncol==0) ) )then
            call debug( (/nrow, ncol, nnz/), "nrow ncol, nnz =", tag=dALLWAYS )
            call stop_program("In set_csr: if nnz is nonzero, nrow and ncol must all be nonzero")
        end if

        if( (nrow==0).eqv.(ncol==0) )then !all zero or nonzero
            if( (nrow==0).and.(ncol==0) )then ! all zero
                ll_reset = .true.
            else !both nonzero
                if(A%initialized)then
                    !if current size different from required size
                    if( (A%nrow /= nrow).or.(A%ncol /= ncol).or.(A%nnz /= nnz) )then
                        ll_reset = .true.
                        ll_allocate = .true.
                    end if
                else !not initialized
                    ll_allocate = .true.
                end if
            end if
        else
            call debug( (/nrow, ncol/), "nrow ncol =", tag=dALLWAYS )
            call stop_program("In set_csr: nrow ncol should all be zero or nonzero")
        end if

        if(ll_reset)call reset_csr(A) !free space if necessary

        if( ll_allocate )then
            A%nrow = nrow
            A%ncol = ncol
            A%nnz  = nnz
            allocate( A%val(A%nnz), A%col_ind(A%nnz), A%row_ptr(A%nrow+1))
        end if

        if( nnz>0 )then !all nonzero
            A%val = 0.0_cp
            A%col_ind = 0
            A%row_ptr = A%nnz+1
            A%initialized = .true.
        end if

        !call debug_state( td_state, tag=dALLWAYS )
    end subroutine set_csr

    !> @brief Reset crs data structure, free dynamically allocated memory
    !! @param [in, out] A sparse matrix data structure
    !<
    subroutine reset_csr(A)
        type(CSR_matrix), intent(inout) :: A

        if(A%initialized) then
            deallocate( A%val, A%col_ind, A%row_ptr )
        end if
        nullify(A%val, A%col_ind, A%row_ptr)

        A%nrow = 0
        A%ncol = 0
        A%nnz  = 0

        A%initialized = .false.
    end subroutine reset_csr

    !> @brief transpose a CSR matrix
    !! @param [in] A sparse matrix data structure to transpose
    !! @param [in] B sparse matrix data structure for the transposed matrix
    !! @details
    !!   this is not an in-place operation.
    !!   Developper notes
    !!   this subroutine exploits some properties of operations on CSR
    !!   matrix to limit the number of variables. For example, the row_ptr
    !!   of the transposed matrix is used for three different purposes:
    !!     - count the number of nonzero of each row of the transposed
    !!     - mark the position in the data array of the next avaible position
    !!        for each row of the transpose
    !!     - its actual role of row_ptr
    !! Something to keep in mind (reason why I implemented this function):
    !! in fortran (pgf90 14.6) the code
    !! y = matmul( transpose(A), x)
    !! and the code
    !! B = transpose(A); y = matmul(B,x)
    !! do not produce the same results.
    !<
    subroutine csr_transpose(A, B)
        type(CSR_matrix), intent(in) :: A
        type(CSR_matrix), intent(inout) :: B
        !local variables
        integer :: i, j, k, q, l

        call new_csr(B, A%ncol, A%nrow, A%nnz)

        ! Step 1
        ! counting the number of nonzeros in each row of the transpose
        B%row_ptr(1:B%nrow) = 0
        do i = 1, A%nrow
            do k = A%row_ptr(i), A%row_ptr(i+1)-1
                j = A%col_ind(k)
                B%row_ptr(j) = B%row_ptr(j) + 1
            end do
        end do
        ! At this point, it is mandatory that sum (B%row_ptr(1:B%nrow))=B%nnz
        if( sum(B%row_ptr(1:B%nrow))/=B%nnz )call stop_program("In csr_transpose, step 1 went wrong")

        ! Step 2
        ! Compute the location in the data vector of the first element
        ! of each row of the transpose
        do i = B%nrow, 1, -1
            B%row_ptr(i) = B%row_ptr(i+1) - B%row_ptr(i)
        end do
        ! At this point, it is mandatory that B%row_ptr(1)=1
        if(B%row_ptr(1)/=1)call stop_program("In csr_transpose, step 2 went wrong")

        ! Step 3
        ! Actual copy
        do i = 1,A%nrow
            do k = A%row_ptr(i), A%row_ptr(i+1)-1
                q = A%col_ind(k)
                l = B%row_ptr(q) ! next index for the q row
                B%col_ind(l) = i ! colums of B are rows of A
                B%val(l) = A%val(k)
                B%row_ptr(q) = B%row_ptr(q) + 1
            end do
        end do
        ! At this point B%row_ptr contains pointer to the begining of the next row

        ! Step 4, re-shift B%row_ptr to actual row pointers
        do i = B%nrow, 1, -1
            B%row_ptr(i+1) = B%row_ptr(i)
        end do
        B%row_ptr(1)=1
        ! At this point, it is mandatory that B%row_ptr(B%nrow+1)=B%nnz+1
        if(B%row_ptr(B%nrow+1)/=B%nnz+1)call stop_program("In csr_transpose, step 4 went wrong")
    end subroutine csr_transpose

    !> @brief multiply a CSR matrix A by a vector.
    !! @param[in] A matrix under the CSR form
    !! @param[in] x input vector, must be of size A%ncol
    !! @param[out] y outut vector, must be of size A%nrow
    !! @param[in] trans logical transpose flag
    !! @details
    !!   There is no check of the dimensions
    !!   if trans flag is present and set to .true. then
    !!   y = A^t x, the transpose is used
    !!   otherwise
    !!   y = Ax
    !<
    subroutine CSR_matvec(A, x, y, trans)
        type(CSR_matrix), intent(in) :: A
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(in out) :: y
        logical, optional :: trans
        !local variables
        integer :: i, k
        real(dp) :: t
        logical :: ll_trans

        if( present(trans) )then
            ll_trans = trans
        else
            ll_trans = .false.
        end if

        if(ll_trans) then
            !apply the transposed matrix
            y = 0.0_dp
            do i = 1, A%nrow
                do k = A%row_ptr(i), A%row_ptr(i+1)-1
                    y( A%col_ind(k) ) = y( A%col_ind(k) ) + x( i )*A%val(k)
                end do
            end do
        else
            do i = 1, A%nrow
                !
                !  Compute the inner product of row i with vector X.
                !
                t = 0.0_dp
                do k = A%row_ptr(i), A%row_ptr(i+1)-1
                    t = t + A%val(k) * x( A%col_ind(k) )
                end do
                y(i) = t

            end do
        end if
    end subroutine CSR_matvec

end module sparse_matrix